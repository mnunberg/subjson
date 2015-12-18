
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

// Includes the source code for the parser itself
#define INCLUDE_JSONSL_SRC
#include "jsonsl_header.h"
#include "subdoc-api.h"
#include "multiget.h"
#include "hkesc.h"
#include <vector>

using namespace Subdoc;

// Utility function to get the number of children
static size_t
get_num_children(const jsonsl_state_st *st)
{
    if (st->type == JSONSL_T_LIST) {
        return st->nelem;
    } else if (st->type == JSONSL_T_OBJECT) {
        return st->nelem / 2;
    } else {
        return 0;
    }
}

namespace {
struct Context;

/** Runtime information for a given path-state pair */
class StateInfo {
private:
    // Parent context object
    Context *m_parent;

    // Path (for matching)
    const Path *m_path;

    // Result (match info is here)
    BaseMatch *m_result;

    // Index (in order to declare the match as 'done')
    size_t m_index;

    std::vector<int16_t> m_levelinfo;

    inline void set_complete();
    inline BaseMatch &result() { return *m_result; }

public:
    StateInfo(Context *cx, size_t index, const Path *path, BaseMatch *match)
    : m_parent(cx), m_path(path), m_result(match), m_index(index) {

        // For the root
        m_levelinfo.resize(Subdoc::Limits::PARSER_DEPTH);
        m_levelinfo[0] = JSONSL_MATCH_POSSIBLE;
        for (size_t ii = 1; ii < m_levelinfo.size(); ii++) {
            // Turn off matching for all subsequent levels until
            // we know otherwise
            m_levelinfo[ii] = JSONSL_MATCH_NOMATCH;
        }
    }

    inline void on_push(jsonsl_state_st *, jsonsl_state_st *);
    inline void on_pop(jsonsl_state_st *, jsonsl_state_st *);
};

struct Context {
    const char *origbuf = NULL;
    jsonsl_t jsn = NULL;
    jsonsl_error_t rc = JSONSL_ERROR_SUCCESS;

    // How many StateInfo objects are still valid?
    size_t num_remaining;

    // Quick lookup of states which are done
    std::vector<int> done_map;

    std::vector<StateInfo> infos;
    HashKey hk;

    bool is_pending(size_t ix) const {
        return !done_map[ix];
    }

    Context(const char *buf, const MultiPath& paths, BaseMatch *results)
    : origbuf(buf) {

        size_t num = paths.size();
        done_map.resize(paths.size());
        num_remaining = paths.size();

        for (size_t ii = 0; ii < paths.size(); ++ii) {
            if (!paths.error_for(ii).success()) {
                done_map[ii] = true;
                num_remaining--;
            }
        }

        // So the pointer locations don't move
        for (size_t ii = 0; ii < num; ++ii) {
            infos.push_back(StateInfo(this, ii, &paths[ii], results + ii));
        }
    }

    static Context *get(jsonsl_t jsn) {
        return reinterpret_cast<Context*>(jsn->data);
    }
};

void
StateInfo::on_push(jsonsl_state_st *parent, jsonsl_state_st *state)
{
    BaseMatch *m = &result();
    if (parent != NULL && m_levelinfo[parent->level] == JSONSL_MATCH_NOMATCH) {
        m_levelinfo[state->level] = JSONSL_MATCH_NOMATCH;
        return;
    }

    if (m->matchres == JSONSL_MATCH_COMPLETE && state->level != m->match_level) {
        // TODO: automatically exclude ourselves for callbacks?
        return; // A child match, perhaps for another path?
    }

    const char *hk;
    size_t nk;
    hk = m_parent->hk.get_hk(nk);

    auto *jpr_p =
            const_cast<jsonsl_jpr_st*>(static_cast<const jsonsl_jpr_st*>(m_path));
    m->matchres = jsonsl_path_match(jpr_p, parent, state, hk, nk);
    m_levelinfo[state->level] = m->matchres;

    if (m->matchres == JSONSL_MATCH_POSSIBLE || m->matchres == JSONSL_MATCH_COMPLETE) {
        m->match_level = state->level;
    } else if (m->matchres == JSONSL_MATCH_NOMATCH) {
        // Indicate that this current state should not be called further for
        // the given tree..
    } else if (m->matchres == JSONSL_MATCH_TYPE_MISMATCH) {
        // Mark as complete, but use an error path
        set_complete();
    } else {
        throw std::runtime_error("Unknown match status received");
    }
}

void
StateInfo::on_pop(jsonsl_state_st *parent, jsonsl_state_st *state)
{
    if (m_levelinfo[state->level] != JSONSL_MATCH_COMPLETE &&
            m_levelinfo[state->level] != JSONSL_MATCH_POSSIBLE) {
        // Not completed, but this is not the a tree we care about
        return;
    }

    BaseMatch *m = &result();

    if (state->level != m->match_level) {
        fprintf(stderr, "CUR: %d. MATCH=%d. INFO=%d\n", state->level, m->match_level, m_levelinfo[state->level]);
        throw std::runtime_error("Should only be called with POSSIBLE/COMPLETE and last state");
    }

    // Only care about marking the results if the current state is
    // indeed the deepest most match (either possible or complete)
    set_complete();
    size_t pos_cur = m_parent->jsn->pos;

    m->type = state->type;
    m->sflags = state->special_flags;
    m->loc_deepest.assign(
        m_parent->origbuf + state->pos_begin,
        pos_cur - state->pos_begin);

    if (state->type == JSONSL_T_STRING || JSONSL_STATE_IS_CONTAINER(state)) {
        // Include trailing tokens, if any
        m->loc_deepest.length++;
    }

    m->num_children = get_num_children(state);
    if (parent != NULL) {
        m->num_siblings = get_num_children(parent);
        if (m->matchres == JSONSL_MATCH_COMPLETE) {
            m->num_siblings--;
            m->position = m->num_siblings;
        }
    }

    if (m->matchres == JSONSL_MATCH_POSSIBLE) {
        if (m->match_level == m_path->size()-1) {
            m->immediate_parent_found = 1;
        }
    } else {
        m->immediate_parent_found = 1;
    }
}

void
StateInfo::set_complete()
{
    m_parent->done_map[m_index] = true;
    if (! --m_parent->num_remaining) {
        jsonsl_stop(m_parent->jsn);
    }
}

static void
push_callback(jsonsl_t jsn, jsonsl_action_t, jsonsl_state_st *state,
    const jsonsl_char_t *at) {
    Context *cx = Context::get(jsn);

    // For push, we ignore any HKEY
    if (state->type == JSONSL_T_HKEY) {
        cx->hk.set_hk_begin(state, at);
        return;
    }

    auto parent = jsonsl_last_state(jsn, state);
    for (size_t ii = 0; ii < cx->infos.size(); ++ii) {
        if (cx->is_pending(ii)) {
            cx->infos[ii].on_push(parent, state);
        }
    }
}

static void
pop_callback(jsonsl_t jsn, jsonsl_action_t, jsonsl_state_st *state,
        const jsonsl_char_t *) {
    Context *cx = Context::get(jsn);

    if (state->type == JSONSL_T_HKEY) {
        cx->hk.set_hk_end(state);
        return;
    }

    auto parent = jsonsl_last_state(jsn, state);
    for (size_t ii = 0; ii < cx->infos.size(); ++ii) {
        if (cx->is_pending(ii)) {
            cx->infos[ii].on_pop(parent, state);
        }
    }
}

static int
error_callback(jsonsl_t jsn, jsonsl_error_t err, jsonsl_state_st*, jsonsl_char_t *)
{
    Context::get(jsn)->rc = err;
    return 0;
}
} // anon namespace

int
MultiMatch::exec_matches(const Loc& loc, const MultiPath& paths, BaseMatch *results, jsonsl_t jsn)
{
    jsonsl_enable_all_callbacks(jsn);
    jsn->action_callback_POP = pop_callback;
    jsn->action_callback_PUSH = push_callback;
    jsn->error_callback = error_callback;

    Context cx(loc.at, paths, results);
    jsn->data = &cx;
    cx.jsn = jsn;

    jsonsl_feed(jsn, loc.at, loc.length);
    return cx.rc;
}

int
MultiMatch::exec_matches(const Loc &doc, const MultiPath& paths, BaseMatch *results)
{
    jsonsl_t jsn = jsonsl_new(Subdoc::Limits::PARSER_DEPTH);
    int rv = exec_matches(doc, paths, results, jsn);
    jsonsl_destroy(jsn);
    return rv;
}

bool
MultiPath::has_overlaps() const
{
    for (size_t ii = 0; ii < m_paths.size(); ++ii) {
        for (size_t jj = 0; jj < m_paths.size(); ++jj) {
            if (jj == ii) {
                // Same path!
                continue;
            }
            if (m_paths[ii].is_related(m_paths[jj])) {
                return true;
            }
        }
    }
    return false;
}
