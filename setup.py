try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, find_packages


import os, sys


NAME = 'deminf_data'

SUPPORTED_PYTHON_VERSIONS = ['3.6', '3.7', '3.8']


# Check python version
if sys.version[0:3] not in SUPPORTED_PYTHON_VERSIONS:
    sys.stderr.write("Python version " + sys.version[0:3] + " is not supported!\n" +
          "Supported versions are " + ", ".join(SUPPORTED_PYTHON_VERSIONS) + "\n")
    sys.stderr.flush()
    sys.exit(1)


# Load up the description from README.md
with open('README.md') as f:
    DESCRIPTION = f.read()

requirements = ['numpy']

# save all data files
all_data_files = list()
for dirname in os.listdir("."):
    if dirname[0].isdigit() and len(dirname.split("_")) >= 4:
        filenames = []
        for filename in ["demographic_model.py",
                         "fs_data.fs",
                         "main_script.py"]:
            filepath = os.path.join(dirname, filename)
            if os.path.exists(filepath):
                filenames.append(filepath)
        if len(filenames) > 0:
            all_data_files.append((dirname, filenames))

setup(
    name=NAME,
    author='Ekaterina Noskova',
    author_email='ekaterina.e.noskova@gmail.com',
    url='https://github.com/noscode/demographic_inference_data',
    description='Data for demographic inference',
    long_description=DESCRIPTION,
    long_description_content_type='text/markdown',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development',
    ],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'deminf_data': ['*.py']
    },
    data_files=all_data_files,
    install_requires=requirements,
)
