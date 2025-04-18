# (optional) To make the project pip-installable and ensure Conda can manage Python dependencies, create a setup.py file in the root directory:
from setuptools import setup, find_packages

setup(
    name='icecream_pkg',
    version='2.0.0',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'icecream_pkg': [
            'resource/*',
            'resource/conjugation/*',
            'resource/oriT_seq/*',
            'example/*',
        ],
    },
    entry_points={
        'console_scripts': [
            'ICEcream=icecream_pkg.ICEcream:main' ,
            'draw_segment=icecream_pkg.plotting_script:plotting_script',
            'arrange_orf=icecream_pkg.organize_orfs:organize_orfs',
            'gbff2gbk=icecream_pkg.gbff2gbk:main' ,
        ],
    },
    install_requires=[
        'biopython',  
    ],
    author='Xiaole (Charlotte) Yin',
    author_email='yinlele99@gmail.com',
    description='ICE identification, classification, visulization',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/your_project',  # Replace with your project's URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Replace with your license
        'Operating System :: OS Independent',
    ],
    python_requires='==3.10.14',  # Specify your Python version
)
