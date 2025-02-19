from setuptools import setup, find_packages

setup(
    name="mrs_pipe",  
    author="Diego Ramírez González", 
    author_email="diegogz95@gmail.com",  
    description="nipype interfaces for MRS software, and workflows for preprocessing svs.",
    # long_description=open("README.md").read(),  
    # long_description_content_type="text/markdown",  
    packages=find_packages(), 
    install_requires=[
        "fsl_mrs>=2.1.13",  
        "nipype>=1.8.0",
    ],
    python_requires=">=3.11"
)
