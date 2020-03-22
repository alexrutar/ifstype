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
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'IfsType'
copyright = '2019, Alex Rutar'
author = 'Alex Rutar'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'sphinx.ext.autodoc',
        #  'sphinx_autodoc_typehints',
        'sphinx.ext.viewcode',
        #  'sphinx.ext.githubpages'
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

source_suffix = '.rst'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = 'alabaster'
html_theme_options = {
    #  'logo': 'gj-logo.png',
    #  'logo_name': True,
    #  'logo_text_align': 'center',
    #  'analytics_id': 'UA-19364636-2',
    #  'show_powered_by': False,
    'show_related': True,
    'github_user': 'alexrutar',
    'github_repo': 'IfsType'
}
html_sidebars = {
    '**': [
        'about.html',
        'localtoc.html',
        "navigation.html",
        'relations.html',
        'searchbox.html',
    ]
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
