import os
import subprocess
import setuptools
from setuptools.command.install import install
from setuptools.command.develop import develop


class SetupPackage(object):
    def check_r(self):
        print("Detecting R installation (and required packages).")
        p = subprocess.Popen([
            "Rscript",
            "-e", "library(rzmq)",
            "-e", "library(rjson)",
        ])

        ret = p.wait()

        if ret != 0:
            raise Exception(
                "R or required packages are not installed. Make sure that "
                "Rscript is in the path and that rjson and rzmq are "
                "installed."
            )

        print("Rscript in path and required packages installed.")

    def make_dealer(self):
        dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "ukbphewas")
        )

        p = subprocess.Popen([
            "make", "-C", dir
        ])

        ret = p.wait()

        if ret != 0:
            raise Exception(
                "Could not compile the ZMQ dealer. Make sure that gcc and "
                "libzmq are installed."
            )

    def run(self):
        self.check_r()
        self.make_dealer()


class CustomInstall(SetupPackage, install):
    def run(self):
        SetupPackage.run(self)
        install.run(self)


class CustomDevelop(SetupPackage, develop):
    def run(self):
        SetupPackage.run(self)
        develop.run(self)


def setup_packages():
    setuptools.setup(
        name="ukbphewas",
        version="2.0.0",
        author="Marc-André Legault",
        author_email="marc-andre.legault.1@umontreal.ca ",
        description="Package to run PheWAS analyses in the UK Biobank.",
        url="https://github.com/legaultmarc",
        packages=setuptools.find_packages(),
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Operating System :: Unix",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        python_requires='>=3.6',
        include_package_data=True,
        install_requires=["pyarrow>=0.15.1", "pandas>=0.25.2",
                          "pyzmq>=18.1.0"],
        zip_safe=False,
        cmdclass={
            "install": CustomInstall,
            "develop": CustomDevelop,
        },
        entry_points={
            "console_scripts": [
                "ukbphewas=ukbphewas.conductor:main"
            ]
        },
    )


if __name__ == "__main__":
    setup_packages()
