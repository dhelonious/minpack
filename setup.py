#!/usr/bin/env python
# encoding: utf-8

from os.path import join


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    dependencies_src = [join("dependencies", "*.f")]

    config = Configuration("minpack", parent_package, top_path)
    config.add_library("dependencies", sources=dependencies_src)
    config.add_extension(
        "_minpack",
        sources=["minpack.pyf", "minpack.f"],
        libraries=["dependencies"],
        include_dirs=["dependencies"],
        depends=dependencies_src
    )

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(
        version="0.0.1",
        description="minpack hybrd and hybrd1 function",
        author="Daniel Dizdarevic",
        author_email="daniel.dizdarevic@gmx.de",
        license="MIT",
        **configuration(top_path="").todict()
    )
