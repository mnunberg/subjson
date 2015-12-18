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

#ifndef SUBDOC_MULTIGET_H
#define SUBDOC_MULTIGET_H
#include "match.h"
#include "loc.h"
#include "path.h"
#include <map>
#include <vector>

namespace Subdoc {

class MultiPath {
public:
    MultiPath() {}

    Error add(const char *s) {
        m_paths.push_back(Path());
        Error rv = Path::rv_to_error(m_paths.back().parse(s));
        if (!rv.success()) {
            m_failed[m_paths.size()-1] = rv;
        }
        return rv;
    }
    Error add(const std::string &s) { return add(s.c_str()); }

    const Path& operator[](size_t idx) const { return m_paths[idx]; }
    size_t size() const { return m_paths.size(); }
    const std::vector<Path>& paths() const { return m_paths; }

    Error error_for(size_t idx) const {
        auto rv = m_failed.find(idx);
        if (rv == m_failed.end()) {
            return Error::SUCCESS;
        }
        return rv->second;
    }

    // Determine if there are overlapping paths. This might be an error for
    // certain mutation operations.
    bool has_overlaps() const;

    void clear() {
        m_failed.clear();
        m_paths.clear();
    }

private:
    std::map<size_t, Error> m_failed;
    std::vector<Path> m_paths;
};

class MultiMatch {
public:
    int exec_matches(const Loc& doc, const MultiPath& paths, BaseMatch *results, jsonsl_t jsn);
    int exec_matches(const Loc& doc, const MultiPath& paths, BaseMatch *results);
};

class LookupOperation {
public:
    void clear();
    LookupOperation();
    ~LookupOperation();

    std::vector<Path> m_paths;
    std::vector<BaseMatch> m_matches;
};

} // namespace

#endif
