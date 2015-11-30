from os import listdir
from os.path import split, isdir, join, isabs
from imp import get_suffixes

import nzmath
_NZMATH_root = split(nzmath.__path__[0])[0]
del nzmath

def import_all(directory):
    """
    Execute 'from (directory) import (all files except "all.py")'.
    """
    try:
        exec('import ' + directory)

        import_files = eval(directory + ".__all__")
        for import_file in import_files:
            if import_file != "all":
                exec('import ' + directory + '.' + import_file)
    except AttributeError: # __all__ does not exist
        pass
    except ImportError: 
        # 1. __init__.py does not exist
        # 2. import_file is not file name
        pass

def search_subdir(directory):
    """
    Return all (1-level) subdirectory-names in directory.
    directory should be a path string.
    """
    return [subdir for subdir in listdir(directory) 
            if isdir(join(directory, subdir))]

def path_to_import_str(path, import_root_path='/'):
    """
    Return '.' separated form (for Python import statements)
    from an (absolute) path.
    """
    for (suffix,mode,imp_type) in get_suffixes():
        dot_idx = path.find(suffix)
        if path.find(suffix) >= 0:
            path = path[:len(suffix)]

    if path.find(import_root_path) >= 0:
        path = path[len(import_root_path):]
    if isabs(path):
        path = path[1:]

    head, tail = split(path)
    if not(head):
        return tail
    head = path_to_import_str(head, import_root_path)
    return head + '.' + tail

def recursive_import_all(rootdir, import_root_path):
    """
    Execute 'from (dir) import *' at all directories in rootdir recursively.
    All directories are searched by using Breadth-First Search 
    from rootdir (rootdir is included in). 

    rootdir should be an absolute path under import_root_path. 
    import_root_path should be a directory 
    whose Python can be searchable on importing process.
    """
    queue = [rootdir]
    while queue:
        target_dir = queue.pop(0)
        import_all(path_to_import_str(target_dir, import_root_path))
        queue.extend([join(target_dir, subdir) 
            for subdir in search_subdir(target_dir)])

script_dir = split(__file__)[0]
recursive_import_all(script_dir, _NZMATH_root)

del listdir, split, isdir, join, isabs, get_suffixes
