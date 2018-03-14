# this new link path policy has fixed ldpath problems
# experienced with out of dyld_load_path libraries on
# OSX without bad effect on other platforms.
cmake_policy(SET CMP0003 NEW)

subdirs(src)
