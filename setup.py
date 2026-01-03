from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        'seq_align',
        ['src/NanoRepeat/seq_align_py.cpp'],
        extra_compile_args=['-O3', '-std=c++11'],
    ),
]

setup(
    packages=find_packages("src"),
    package_dir={"": "src"},
    data_files=[("", ["LICENSE"])],
    scripts=['src/NanoRepeat/nanoRepeat.py', 'src/NanoRepeat/nanoRepeat-joint.py'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)