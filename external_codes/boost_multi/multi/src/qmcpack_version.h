// Wrapper around the auto-generated Git repository revision
// information file (git-rev.h)
// If not building from a git repository, the git-rev.h file is empty

#ifndef QMCPACK_VERSION_INCLUDE
#define QMCPACK_VERSION_INCLUDE

#define STR_EXPAND(x) #x
#define STR(x) STR_EXPAND(x)

#ifdef IS_GIT_PROJECT
#include "git-rev.h"
#endif

#ifdef GIT_BRANCH_RAW
#define QMCPACK_GIT_BRANCH STR(GIT_BRANCH_RAW)
#define QMCPACK_GIT_HASH STR(GIT_HASH_RAW)
#define QMCPACK_GIT_COMMIT_LAST_CHANGED STR(GIT_COMMIT_LAST_CHANGED_RAW)
#define QMCPACK_GIT_COMMIT_SUBJECT GIT_COMMIT_SUBJECT_RAW
#endif

#endif
