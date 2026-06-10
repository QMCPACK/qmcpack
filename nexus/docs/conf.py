# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Nexus'
copyright = '2020, Jaron T. Krogel'
author = 'Jaron T. Krogel'

pygments_style = "default"
highlight_language = "python"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
import sys
from pathlib import Path
sys.path.append(Path(__file__).parent.parent.resolve())

from intersphinx_registry import get_intersphinx_mapping

extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "pydata_sphinx_theme",
    "sphinx_design",
    "sphinx_copybutton",
    "numpydoc",
]

copybutton_exclude = '.linenos, .gp, .go' # Don't copy line numbers, prompts, or outputs.

bibtex_bibfiles = ['bibs/methods.bib']

numpydoc_class_members_toctree = False
numpydoc_show_class_members = True
numpydoc_xref_ignore = {"optional", "type_without_description"}
numpydoc_xref_param_type = True
numpydoc_xref_aliases = {
    "ArrayLike": "numpy.typing.ArrayLike",
    "NDArray": "numpy.typing.NDArray",
    "array_like": "numpy.typing.ArrayLike",
    "ndarray": "numpy.typing.NDArray",
}
numfig = True

add_function_parentheses = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]
html_theme_options = {
    "logo": {
        #"alt_text": "Nexus Docs",
        "text": "Nexus Documentation",
        "image_dark": "_static/nexus_logo.svg",
        "image_light": "_static/nexus_logo.svg",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/QMCPACK/qmcpack",
            "icon": "fa-brands fa-github",
        },
    ],
    "navbar_end": [
        "search-button",
        "theme-switcher",
        "navbar-icon-links",
    ],
    "show_toc_level": 3,
    "collapse_navigation": True,
    "secondary_sidebar_items": ["page-toc"],
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",
}
html_favicon = "_static/nexus_logo.svg"

html_context = {"default_mode": "auto"}

intersphinx_mapping = get_intersphinx_mapping(packages=["python", "numpy"])
