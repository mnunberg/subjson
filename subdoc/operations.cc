/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/*
*     Copyright 2015 Couchbase, Inc
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#define INCLUDE_SUBDOC_STRING_SRC
#define INCLUDE_SUBDOC_NTOHLL
#define NOMINMAX // For Visual Studio

#include "operations.h"
#include "validate.h"
#include "util.h"
#include <errno.h>
#include <inttypes.h>
#include <string>
#include <limits>

using Subdoc::Loc;
using Subdoc::Error;
using Subdoc::Operation;
using Subdoc::Path;
using Subdoc::Match;
using Subdoc::Command;

static Loc loc_COMMA(",", 1);
static Loc loc_QUOTE("\"", 1);
static Loc loc_COMMA_QUOTE(",\"", 2);
static Loc loc_QUOTE_COLON("\":", 2);

/**
 * Performs common matching using the currently designated path. Note that
 * unlike do_get(), a missing match is not an error. This is the reason do_get()
 * and do_match_common() are separate functions, as do_get() (i.e. a simple get
 * operation) considers a non-complete match as a failure.
 *
 * @return success if no parse error or mismatch happened, failure otherwise.
 */
Error
Operation::do_match_common(Match::SearchOptions options)
{
    m_match.extra_options = options;
    m_match.exec_match(m_doc, m_path, m_jsn);

    if (m_match.matchres == JSONSL_MATCH_TYPE_MISMATCH) {
        return Error::PATH_MISMATCH;
    } else if (m_match.status != JSONSL_ERROR_SUCCESS) {
        if (m_match.status == JSONSL_ERROR_LEVELS_EXCEEDED) {
            return Error::DOC_ETOODEEP;
        } else {
            return Error::DOC_NOTJSON;
        }
    } else {
        return Error::SUCCESS;
    }
}

Error
Operation::do_get()
{
    if (m_match.matchres != JSONSL_MATCH_COMPLETE) {
        return Error::PATH_ENOENT;
    }

    m_result->m_match = m_match.loc_deepest;
    return Error::SUCCESS;
}

/* Start at the beginning of the buffer, stripping first comma */
#define STRIP_FIRST_COMMA 1

/* Start at the end of the buffer, stripping last comma */
#define STRIP_LAST_COMMA 2

static void
strip_comma(Loc *loc, int mode)
{
    unsigned ii;
    if (mode == STRIP_FIRST_COMMA) {
        for (ii = 0; ii < loc->length; ii++) {
            if (loc->at[ii] == ',') {
                loc->at += (ii + 1);
                loc->length -= (ii + 1);
                return;
            }
        }
    } else {
        for (ii = loc->length; ii; ii--) {
            if (loc->at[ii-1] == ',') {
                loc->length = ii-1;
                return;
            }
        }
    }
}

Error
Operation::do_remove()
{
    Error rv = do_match_common(Match::GET_FOLLOWING_SIBLINGS);

    if (!rv.success()) {
        return rv;
    } else if (m_match.matchres != JSONSL_MATCH_COMPLETE) {
        return Error::PATH_ENOENT;
    }

    // Deepest is the actual match
    const Loc& match_loc = m_match.loc_deepest;
    /*
     * If match is the _last_ item, we need to strip the _last_ comma from
     * doc[0] (since doc[1] will never have a trailing comma); e.g.
     *  NEWDOC[0] = { ..., ..., ... (---> , <---)
     *  [deleted]
     *  NEWDOC[1] = } // Never any commas here!
     *
     * If match is the _first_ item, we need to strip the _first_ comma
     * from newdoc[1]:
     *
     *  NEWDOC[0] = { // never any commas here!
     *  [deleted]
     *  NEWDOC[1] = (---> , <---) ... , ... }
     *
     *
     *
     * If match is any other item, we need to strip the _last_ comma from
     * doc[0], which is the comma that preceded this item.
     */

    /* Remove the matches, starting from the beginning of the key */
    if (!m_match.loc_key.empty()) {
        newdoc_at(0).end_at_begin(m_doc, m_match.loc_key, Loc::NO_OVERLAP);
    } else {
        newdoc_at(0).end_at_begin(m_doc, match_loc, Loc::NO_OVERLAP);
    }

    newdoc_at(1).begin_at_end(m_doc, match_loc, Loc::NO_OVERLAP);

    if (m_match.num_siblings) {
        if (m_match.is_last()) {
            /*
             * NEWDOC[0] = [a,b,c, <-- Strip here
             * MATCH     = d
             * NEWDOC[1] = ]
             */
            strip_comma(&newdoc_at(0), STRIP_LAST_COMMA);
        } else {
            /*
             * NEWDOC[0] = [a,b,
             * MATCH     = c
             * NEWDOC[1] = Strip here -->, d]
             */
            strip_comma(&newdoc_at(1), STRIP_FIRST_COMMA);
        }
    }

    m_result->m_newlen = 2;
    return Error::SUCCESS;
}

Error
Operation::do_store_dict()
{
    // Check that the last element is not an array first.
    if (m_optype != Command::REPLACE &&
            path()[path().size()-1].ptype == JSONSL_PATH_NUMERIC) {
        return Error::PATH_EINVAL;
    }
    if (m_match.matchres != JSONSL_MATCH_COMPLETE) {
        if (m_optype == Command::REPLACE) {
            // Nothing to replace
            return Error::PATH_ENOENT;
        }
        if (!m_match.immediate_parent_found) {
            if (m_optype.is_mkdir_p()) {
                return do_mkdir_p(MKDIR_P_DICT);
            } else {
                return Error::PATH_ENOENT;
            }
        }
        // Otherwise, the immediate parent is found, and we can
        // continue to add/upsert
    } else if (m_match.matchres == JSONSL_MATCH_COMPLETE) {
        if (m_optype.base() == Command::DICT_ADD) {
            return Error::DOC_EEXISTS;
        }
    }

    if (m_match.matchres == JSONSL_MATCH_COMPLETE) {
        const Loc& match_loc = m_match.loc_deepest;
        /* 1. Remove the old value from the first segment */
        newdoc_at(0).end_at_begin(m_doc, match_loc, Loc::NO_OVERLAP);

        /* 2. Insert the new value */
        newdoc_at(1) = m_userval;

        /* 3. Insert the rest of the document */
        newdoc_at(2).begin_at_end(m_doc, match_loc, Loc::NO_OVERLAP);
        m_result->m_newlen = 3;

    } else if (m_match.immediate_parent_found) {
        // Deepest is parent
        const Loc& parent_loc = m_match.loc_deepest;

        newdoc_at(0).end_at_end(m_doc, parent_loc, Loc::NO_OVERLAP);
        /*TODO: The key might have a literal '"' in it, which has been escaped? */
        if (m_match.num_siblings) {
            newdoc_at(1) = loc_COMMA_QUOTE; /* ," */
        } else {
            newdoc_at(1) = loc_QUOTE /* " */;
        }

        /* Create the actual key: */
        auto& comp = m_path->back();
        newdoc_at(2).at = comp.pstr;
        newdoc_at(2).length = comp.len;

        /* Close the quote and add the dictionary key */
        newdoc_at(3) = loc_QUOTE_COLON; /* ": */
        /* new value */
        newdoc_at(4) = m_userval;
        /* Closing tokens */
        newdoc_at(5).begin_at_end(m_doc, parent_loc, Loc::OVERLAP);
        m_result->m_newlen = 6;
    }
    return Error::SUCCESS;
}

/**
 * This function will inspect the match to find the deepest most parent
 * and add missing entries in the document from the path. It will then
 * insert the given value (either as a value to the dictionary key, or
 * as an array with a single element (being the value) as the value for
 * the dictionary key.
 *
 * In all cases, intermediate non-existing path elements must all refer
 * to dictionary keys, not array elements (since array elements are in
 * sequence there is no course of action if the array is empty and the
 * path refers to a middle element).
 *
 * @param mode How to insert the value
 * @return status
 */
Error
Operation::do_mkdir_p(MkdirPMode mode)
{
    unsigned ii;
    const Loc& parent_loc = m_match.loc_deepest;
    newdoc_at(0).end_at_end(m_doc, parent_loc, Loc::NO_OVERLAP);

    /* doc_new LAYOUT:
     *
     * [0] = HEADER
     * [1] = _P header (bkbuf_extra)
     * [2] = USER TEXT
     * [3] = _P trailer (bkbuf_extra)
     * [4] = TRAILER
     */

    /* figure out the components missing */
    /* THIS IS RESERVED FOR doc_new[1]! */
    if (m_match.num_siblings) {
        m_result->m_bkbuf += ',';
    }

    /* Insert the first item. This is a dictionary key without any object
     * wrapper: */
    const Path::Component* comp = &m_path->get_component(m_match.match_level);
    if (comp->ptype == JSONSL_PATH_NUMERIC) {
        // If it's *not* a dictionary key, don't insert it!
        return Error::PATH_ENOENT;
    }
    m_result->m_bkbuf += '"';
    m_result->m_bkbuf.append(comp->pstr, comp->len);
    m_result->m_bkbuf += "\":";

    /* The next set of components must be created as entries within the
     * newly created key */
    for (ii = m_match.match_level + 1; ii < m_path->size(); ii++) {
        comp = &m_path->get_component(ii);
        if (comp->ptype != JSONSL_PATH_STRING) {
            return Error::PATH_ENOENT;
        }
        m_result->m_bkbuf += "{\"";
        m_result->m_bkbuf.append(comp->pstr, comp->len);
        m_result->m_bkbuf += "\":";
    }
    if (mode == MKDIR_P_ARRAY) {
        m_result->m_bkbuf += '[';
    }

    newdoc_at(1).length = m_result->m_bkbuf.size();

    if (mode == MKDIR_P_ARRAY) {
        m_result->m_bkbuf += ']';
    }

    for (ii = m_match.match_level+1; ii < m_path->size(); ii++) {
        m_result->m_bkbuf += '}';
    }
    newdoc_at(3).length = m_result->m_bkbuf.size() - newdoc_at(1).length;

    /* Set the buffers */
    newdoc_at(1).at = m_result->m_bkbuf.data();
    newdoc_at(3).at = m_result->m_bkbuf.data() + newdoc_at(1).length;
    newdoc_at(2) = m_userval;

    newdoc_at(4).begin_at_end(m_doc, parent_loc, Loc::OVERLAP);
    m_result->m_newlen = 5;

    return Error::SUCCESS;
}

/* Inserts a single element into an empty array */
Error
Operation::insert_singleton_element()
{
    const Loc& array_loc = m_match.loc_deepest;
    /* First segment is ... [ */
    newdoc_at(0).end_at_begin(m_doc, array_loc, Loc::OVERLAP);
    /* User: */
    newdoc_at(1) = m_userval;
    /* Last segment is ... ] */
    newdoc_at(2).begin_at_end(m_doc, array_loc, Loc::OVERLAP);

    m_result->m_newlen = 3;
    return Error::SUCCESS;
}

Error
Operation::do_list_prepend()
{
    Error rv;

    // Search for first element, perform the search, then pop it
    if (m_path->add_array_index(0) != JSONSL_ERROR_SUCCESS) {
        return Error::PATH_E2BIG;
    }

    rv = do_match_common(Match::GET_MATCH_ONLY);
    m_path->pop();

    if (!rv.success()) {
        return rv;
    }

    if (m_match.matchres != JSONSL_MATCH_COMPLETE) {
        if (m_match.immediate_parent_found) {
            return insert_singleton_element();
        } else if (m_optype.is_mkdir_p()) {
            return do_mkdir_p(MKDIR_P_ARRAY);
        } else {
            return Error::PATH_ENOENT;
        }
    }

    /* LAYOUT:
     * NEWDOC[0] [....
     * NEWDOC[1] = USER
     * NEWDOC[2] = ,
     * NEWDOC[3] = REST
     */

    const Loc& orig_firstloc = m_match.loc_deepest;
    newdoc_at(0).end_at_begin(m_doc, orig_firstloc, Loc::NO_OVERLAP);
    newdoc_at(1) = m_userval;
    newdoc_at(2) = loc_COMMA;
    newdoc_at(3).begin_at_begin(m_doc, orig_firstloc);

    m_result->m_newlen = 4;
    return Error::SUCCESS;
}

Error
Operation::do_list_append()
{
    Error rv;

    // Find the match itself, no magic needed!
    bool is_unique = m_optype.base() == Command::ARRAY_ADD_UNIQUE;
    if (is_unique) {
        m_match.ensure_unique = m_userval;
    }

    rv = do_match_common(Match::GET_MATCH_ONLY);
    if (!rv.success()) {
        // Parse error (ENOENT excluded)
        return rv;
    }

    if (m_match.matchres != JSONSL_MATCH_COMPLETE) {
        // Not complete. Determine if mkdir_p should be used
        if (m_optype.is_mkdir_p()) {
            return do_mkdir_p(MKDIR_P_ARRAY);
        } else {
            return Error::PATH_ENOENT;
        }
    }

    // All other errors should be handled above
    SUBDOC_ASSERT(m_match.matchres == JSONSL_MATCH_COMPLETE);

    // Incorrect terminal type!
    if (m_match.type != JSONSL_T_LIST) {
        return Error::PATH_MISMATCH;
    }

    if (is_unique && m_match.unique_item_found) {
        // Unique item already exists
        return Error::DOC_EEXISTS;
    }
    const Loc& array_loc = m_match.loc_deepest;

    if (m_match.num_children == 0) {
        // Special case, no children; append and prepend are the same here!
        return insert_singleton_element();
    }

    /* Last element */
    newdoc_at(0).end_at_end(m_doc, array_loc, Loc::NO_OVERLAP);
    /* Insert comma */
    newdoc_at(1) = loc_COMMA;
    /* User */
    newdoc_at(2) = m_userval;
    /* Parent end */
    newdoc_at(3).begin_at_end(m_doc, array_loc, Loc::OVERLAP);

    m_result->m_newlen = 4;
    return Error::SUCCESS;
}

Error
Operation::do_insert()
{
    auto& lastcomp = m_path->get_component(m_path->size()-1);
    if (lastcomp.ptype != JSONSL_PATH_NUMERIC) {
        return Error::PATH_EINVAL;
    }
    if (lastcomp.is_neg) {
        // Negative paths are invalid for insert operations
        return Error::PATH_EINVAL;
    }

    Error status = do_match_common(Match::GET_MATCH_ONLY);
    if (!status.success()) {
        return status;
    }

    if (m_match.matchres == JSONSL_MATCH_COMPLETE) {
        const Loc& match_loc = m_match.loc_deepest;

        /*
         * DOCNEW[0] = ... [
         * DOCNEW[1] = USER
         * DOCNEW[2] = ,
         * DOCNEW[3] = MATCH
         */
        m_result->m_newlen = 4;
        newdoc_at(0).end_at_begin(m_doc, match_loc, Loc::NO_OVERLAP);
        newdoc_at(1) = m_userval;
        newdoc_at(2) = loc_COMMA;
        newdoc_at(3).begin_at_begin(m_doc, match_loc);
        return Error::SUCCESS;

    } else if (m_match.immediate_parent_found) {
        const Loc& array_loc = m_match.loc_deepest;
        if (m_match.num_siblings == 0 && lastcomp.idx == 0) {
            // Equivalent to prepend/push_first
            /*
             * DOCNEW[0] = ... [
             * DOCNEW[1] = USER
             * DOCNEW[2] = ] ...
             */
            m_result->m_newlen = 3;
            newdoc_at(0).end_at_begin(m_doc, array_loc, Loc::OVERLAP);
            newdoc_at(1) = m_userval;
            newdoc_at(2).begin_at_end(m_doc, array_loc, Loc::OVERLAP);
            return Error::SUCCESS;

        } else if (lastcomp.idx == m_match.num_siblings) {
            // Equivalent to append/push_last
            /*
             * (assume DOC = [a,b,c,d,e]
             * DOCNEW[0] = e (last char before ']')
             * DOCNEW[1] = , (since there are items in the list)
             * DOCNEW[2] = USER
             * DOCNEW[3] = ]
             */
            m_result->m_newlen = 4;
            newdoc_at(0).end_at_end(m_doc, array_loc, Loc::NO_OVERLAP);
            newdoc_at(1) = loc_COMMA;
            newdoc_at(2) = m_userval;
            newdoc_at(3).begin_at_end(m_doc, array_loc, Loc::OVERLAP);
            return Error::SUCCESS;

        } else {
            return Error::PATH_ENOENT;
        }
    } else {
        return Error::PATH_ENOENT;
    }
}

static Error
parse_int64(const Loc& orig, int64_t& outval)
{
    const char *cur = orig.at;
    size_t n = orig.length;

    if (!n) {
        return Error::VALUE_EMPTY; // Empty
    }
    if (*cur == '-') {
        cur++;
        if (!--n) {
            return Error::VALUE_EBADNUMBER;
        }
    }
    if (*cur == '0') {
        // Can't start with a zero. It's either a leading zero or an empty
        // delta
        return Error::VALUE_EZERODELTA;
    }

    outval = 0;
    for (size_t ii = 0; ii < n; ii++) {
        if (!isdigit(cur[ii])) {
            return Error::VALUE_EBADNUMBER; // Not a number!
        }
        // Get the numeric value of the digit, with '0' being the lowest
        // value digit character in the ascii table, and with digits appearing
        // in order, such that '9' (0x39) - '0' (0x30) == 0

        uint64_t newval = (outval * 10) + (cur[ii] - '0');
        if (newval < static_cast<uint64_t>(outval) ||
                newval > std::numeric_limits<int64_t>::max()) {
            // mismatch
            return Error::DELTA_E2BIG;
        }
        outval = static_cast<int64_t>(newval);
    }
    if (*orig.at == '-') {
        outval *= -1;
    }
    return Error::SUCCESS;
}

Error
Operation::do_arith_op()
{
    Error status;
    int64_t delta;
    int64_t numres = 0;
    // Verify the digit first

    status = parse_int64(m_userval, delta);
    if (!status.success()) {
        return status;
    }

    /* Find the number first */
    status = do_match_common(m_optype.is_mkdir_p() ?
            Match::GET_FOLLOWING_SIBLINGS : Match::GET_MATCH_ONLY);

    if (status != Error::SUCCESS) {
        return status;
    }

    if (m_match.matchres == JSONSL_MATCH_COMPLETE) {
        const Loc& num_loc = m_match.loc_deepest;
        if (m_match.type != JSONSL_T_SPECIAL) {
            return Error::PATH_MISMATCH;
        } else if (m_match.sflags & ~(JSONSL_SPECIALf_NUMERIC)) {
            return Error::PATH_MISMATCH;
        } else  {
            errno = 0;
            numres = strtoll(num_loc.at, NULL, 10);

            if (errno == ERANGE) {
                return Error::NUM_E2BIG;
            }

            /* Calculate what to place inside the buffer. We want to be gentle here
             * and not force 64 bit C arithmetic to confuse users, so use proper
             * integer overflow/underflow with a 64 (or rather, 63) bit limit. */
            if (delta >= 0 && numres >= 0) {
                if (std::numeric_limits<int64_t>::max() - delta < numres) {
                    return Error::DELTA_E2BIG;
                }
            } else if (delta < 0 && numres < 0) {
                if (delta < std::numeric_limits<int64_t>::min() - numres) {
                    return Error::DELTA_E2BIG;
                }
            }

            numres += delta;
            m_result->m_numbuf = std::to_string(numres);
        }
    } else {
        if (!m_optype.is_mkdir_p() && !m_match.immediate_parent_found) {
            return Error::PATH_ENOENT;
        }

        if (m_match.type != JSONSL_T_OBJECT) {
            return Error::PATH_ENOENT;
        }

        m_result->m_numbuf = std::to_string(delta);
        m_userval.at = m_result->m_numbuf.data();
        m_userval.length = m_result->m_numbuf.size();
        m_optype = Command::DICT_ADD_P;
        if ((status = do_store_dict()) != Error::SUCCESS) {
            return status;
        }
        m_result->m_match = m_match.loc_deepest = m_userval;
        return Error::SUCCESS;
    }

    // Unlike other ops, we modify the match value itself, as this is an
    // in/out parameter

    /* Preamble */
    newdoc_at(0).end_at_begin(m_doc, m_match.loc_deepest, Loc::NO_OVERLAP);

    /* New number */
    newdoc_at(1).at = m_result->m_numbuf.data();
    newdoc_at(1).length = m_result->m_numbuf.size();

    /* Postamble */
    newdoc_at(2).begin_at_end(m_doc, m_match.loc_deepest, Loc::NO_OVERLAP);
    m_result->m_newlen = 3;

    m_match.loc_deepest.at = m_result->m_numbuf.data();
    m_match.loc_deepest.length = m_result->m_numbuf.size();
    m_result->m_match = m_match.loc_deepest;
    return Error::SUCCESS;
}

Error
Operation::validate(int mode, int depth)
{
    if (!m_userval.empty()) {
        int rv = Validator::validate(m_userval, m_jsn, depth, mode);
        switch (rv) {
        case JSONSL_ERROR_SUCCESS:
            return Error::SUCCESS;
        case JSONSL_ERROR_LEVELS_EXCEEDED:
            return Error::VALUE_ETOODEEP;
        default:
            return Error::VALUE_CANTINSERT;
        }
    } else {
        return Error::VALUE_EMPTY;
    }
}

int
Operation::get_maxdepth(DepthMode mode) const
{
    if (mode == DepthMode::PATH_HAS_NEWKEY) {
        return (Limits::MAX_COMPONENTS + 1) - m_path->size();
    } else {
        return Limits::MAX_COMPONENTS - m_path->size();
    }
}

Error
Operation::op_exec(const char *pth, size_t npth)
{
    int rv = m_path->parse(pth, npth);
    Error status;

    if (rv != 0) {
        if (rv == JSONSL_ERROR_LEVELS_EXCEEDED) {
            return Error::PATH_E2BIG;
        } else {
            return Error::PATH_EINVAL;
        }
    }

    switch (m_optype) {
    case Command::GET:
    case Command::EXISTS:
        status = do_match_common(Match::GET_MATCH_ONLY);
        if (status != Error::SUCCESS) {
            return status;
        }
        return do_get();

    case Command::DICT_ADD:
    case Command::DICT_ADD_P:
    case Command::DICT_UPSERT:
    case Command::DICT_UPSERT_P:
    case Command::REPLACE:
        if (m_path->size() == 1) {
            /* Can't perform these operations on the root element since they
             * will invalidate the JSON or are otherwise meaningless. */
            return Error::VALUE_CANTINSERT;
        }

        status = validate(Validator::PARENT_DICT, get_maxdepth(PATH_HAS_NEWKEY));
        if (!status.success()) {
            return status;
        }
        status = do_match_common(Match::GET_MATCH_ONLY);
        if (!status.success()) {
            return status;
        }
        return do_store_dict();

    case Command::REMOVE:
        if (m_path->size() == 1) {
            // Can't remove root element!
            return Error::VALUE_CANTINSERT;
        }
        return do_remove();

    case Command::ARRAY_PREPEND:
    case Command::ARRAY_PREPEND_P:
        status = validate(Validator::PARENT_ARRAY, get_maxdepth(PATH_IS_PARENT));
        if (!status.success()) {
            return status;
        }
        return do_list_prepend();

    case Command::ARRAY_APPEND:
    case Command::ARRAY_APPEND_P:
    case Command::ARRAY_ADD_UNIQUE:
    case Command::ARRAY_ADD_UNIQUE_P: {
        int validmode = Validator::PARENT_ARRAY;
        if (m_optype.base() == Command::ARRAY_ADD_UNIQUE) {
            // Uniqueness must contain a single, primitive value.
            validmode |= Validator::VALUE_PRIMITIVE | Validator::VALUE_SINGLE;
        }

        status = validate(validmode, get_maxdepth(PATH_IS_PARENT));
        if (!status.success()) {
            return status;
        }
        return do_list_append();
    }

    case Command::ARRAY_INSERT:
        status = validate(Validator::PARENT_ARRAY, get_maxdepth(PATH_HAS_NEWKEY));
        if (!status.success()) {
            return status;
        }
        return do_insert();

    case Command::COUNTER:
    case Command::COUNTER_P:
        // no need to check for depth here, since if the path itself is too
        // big, it will fail during path parse-time
        return do_arith_op();

    default:
        return Error::GLOBAL_ENOSUPPORT;

    }
}

Operation::Operation()
: m_path(new Path()),
  m_jsn(Match::jsn_alloc()),
  m_optype(Command::GET),
  m_result(NULL)
{
}

Operation::~Operation()
{
    clear();
    delete m_path;
    Match::jsn_free(m_jsn);
}

void
Operation::clear()
{
    m_path->clear();
    m_match.clear();
    m_userval.length = 0;
    m_userval.at = NULL;
    m_result = NULL;
    m_optype = Command::GET;
}

Operation *
subdoc_op_alloc()
{
    return new Operation();
}

void
subdoc_op_free(Operation *op)
{
    delete op;
}

/* Misc */
const char *
Error::description() const
{
    switch (m_code) {
    case Error::SUCCESS:
        return "Success";
    case Error::PATH_ENOENT:
        return "Requested path does not exist in document";
    case Error::PATH_MISMATCH:
        return "The path specified treats an existing document entry as the wrong type";
    case Error::PATH_EINVAL:
        return "Path syntax error";
    case Error::DOC_NOTJSON:
        return "The document is not JSON";
    case Error::DOC_EEXISTS:
        return "The requested path already exists";
    case Error::PATH_E2BIG:
        return "The path is too big";
    case Error::NUM_E2BIG:
        return "The number specified by the path is too big";
    case Error::DELTA_E2BIG:
        return "The combination of the existing number and the delta will result in an underflow or overflow";
    case Error::VALUE_CANTINSERT:
        return "The new value cannot be inserted in the context of the path, as it would invalidate the JSON";
    case Error::VALUE_EMPTY:
        return "Expected non-empty value for command";
    case Error::VALUE_ETOODEEP:
        return "Adding this value would make the document too deep";
    case Error::GLOBAL_ENOMEM:
        return "Couldn't allocate memory";
    case Error::GLOBAL_ENOSUPPORT:
        return "Operation not implemented";
    case Error::GLOBAL_UNKNOWN_COMMAND:
        return "Unrecognized command code";
    case Error::VALUE_EBADNUMBER:
        return "String is not a JSON number";
    case Error::VALUE_EZERODELTA:
        return "Delta value is zero";
    case Error::DOC_ETOODEEP:
        return "Document is too deep to parse";
    default:
        return "Unknown error code";
    }
}
