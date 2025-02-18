from setuptools import setup, find_packages

setup(
    name="mrs_pipe",  
    author="Diego Ramírez González", 
    author_email="diegogz95@gmail.com",  
    description="Una breve descripción de tu paquete",
    # long_description=open("README.md").read(),  
    # long_description_content_type="text/markdown",  
    packages=find_packages(), 
    install_requires=[
        "fsl_mrs>=2.2.0",  
        "nipype>=1.8.0",
    ],
    python_requires=">=3.11"
)
