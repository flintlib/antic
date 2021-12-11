r"""
Pip-compatible setup instructions for Antic.

Being a C/C++-library with a few dependencies, Antic cannot be sanely
installed with pip from binary builds. It should be installed through a proper
package manager. When this is not possible, it needs to be installed from
source with configure-make-make-install. However, in some scenarios, hiding
this build inside a pip-package is advantageous, e.g., when installing inside a
SageMath source build or when packages have to be installed with a
requirements.txt.
"""
# ####################################################################
#  This file is part of Antic.
#
#        Copyright (C) 2021 Julian Rüth
#
#  Antic is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 2.1 of the License, or (at your option)
#  any later version.
#
#  Antic is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along with
#  Antic. If not, see <https://www.gnu.org/licenses/>.
# ###################################################################
import os
import inspect
import shutil
from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build
from setuptools.command.sdist import sdist
from subprocess import check_call
from contextlib import contextmanager


@contextmanager
def cwd(path):
    r"""
    Change the current working directory to `path` while inside this context.
    """
    pwd = os.getcwd()
    os.chdir(path)
    try:
        yield path
    finally:
        os.chdir(pwd)


class AutotoolsCommand:
    r"""
    Helpers that provide variables like the ones available to automake
    Makefiles.
    """
    @property
    def builddir(self):
        r"""
        Return the path of the build directory, equivalent to @builddir@ in automake.

        This property is only available when a build has been invoked as part of this setup.py run.
        In particular, this is the directory passed with --build-base.
        Note that "setup.py install" does not accept a --build-base so we have
        to make sure that this happens to be build/ relative to the current
        working directory during the install phase since otherwise the install
        step won't be able to find the assets that have been built by the build
        step.
        """
        if "build" not in self.distribution.command_obj:
            raise ValueError("There is no build directory since no build is being performed in this invocation of setup.py")

        return self.distribution.command_obj["build"].build_base

    @property
    def abs_builddir(self):
        r"""
        Return the absolute path of the build directory, equivalent to @abs_builddir@ in automake.

        The limitations of `builddir` apply to this property as well.
        """
        builddir = self.builddir
        if not os.path.isabs(builddir):
            builddir = os.path.join(self.abs_srcdir, builddir)

        return builddir

    @property
    def destdir(self):
        r"""
        Return the installation prefix for this package in site-packages (or the user directory.)

        This is the value that you want to pass to a configure's --prefix flag.
        Note that naturally this value is only available when setup was asked
        to install this package. In particular, this is not available when
        trying to build a wheel.

        As a consequence this value is also not available initially when
        invoking `pip install` since that tries to build a wheel first. You
        might want to invoke pip install with `--no-binary :all:` so that pip
        skips to a regular install where this value is available.
        """
        if "install" not in self.distribution.command_obj:
            raise ValueError("Cannot determine installation prefix in this build which does not install.")
        return os.path.join(self.distribution.command_obj["install"].install_lib, self.distribution.get_name())

    @property
    def abs_srcdir(self):
        r"""
        Return the absolute path of the source directory, i.e., the directory where this setup.py is.

        This is the equivalent to @abs_srcdir@ in automake.
        """
        return os.path.abspath(os.path.dirname(__file__) or ".")

    @property
    def MAKE(self):
        r"""
        Return the name of the make command which might have been overridden by
        the MAKE environment variable.
        """
        return os.getenv("MAKE", "make")


class NoBinaryBuild(build):
    r"""
    Disables building binary wheels for this package.

    Since such binary wheels are not relocatable because we hard code library
    paths into our hedaer files so cppyy known where things are located.
    """

    def run(self):
        raise NotImplementedError("No binary wheels can be built for Antic currently because the installation prefix is hard-coded in some of its header files. To skip this step when using pip, run with --no-binary :all:")


class ConfigureMakeInstall(install, AutotoolsCommand):
    r"""
    Builds and installs Antic by running configure and make install.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # We do not perform a separate build since we don't know the prefix during such a build yet.
        self.skip_build = True

    def run(self):
        super().run()

        # Perform a build and install into a directory such as
        # site-packages/antic.
        with cwd(self.abs_srcdir):
            check_call([os.path.join(self.abs_srcdir, "configure"), f"--prefix={self.destdir}"])
            check_call([self.MAKE, "install"])

    def get_outputs(self):
        r"""
        Return the installed files/directories so we know how to uninstall Antic again.
        """
        return super().get_outputs() + [self.destdir]


setup(
    name='antic',
    author='the Antic authors',
    url='https://github.com/wbhart/antic',
    # We cannot encode all our dependencies since they are mostly not available
    # as Python packages. Also we want to build with system packages that are
    # detected by the configure script.
    install_requires=[],
    long_description=inspect.cleandoc(r"""
        Antic is an algebraic number theory library in C.

        We do not recommend to install Antic from PyPI as it has dependencies that are not available on PyPI.

        Please consult Antic's home page for further details: https://github.com/wbhart/antic
        """),
    version='0.2.1',
    license='LGPL 2.1+',
    license_files=('LICENSE', 'gpl-2.0.txt'),
    setup_requires=["wheel"],
    cmdclass={
        'build': NoBinaryBuild,
        'install': ConfigureMakeInstall,
    },
)
