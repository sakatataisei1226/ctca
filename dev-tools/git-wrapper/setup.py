import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="git-wrapper", # Replace with your own username
    version="0.1.0",
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'git = wrapper.main:main',
        ],
    },
    author="Nkzono99",
    author_email="",
    description="Git command wrapper",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)