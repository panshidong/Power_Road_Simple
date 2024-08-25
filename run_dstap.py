import jpype
import jpype.imports
from jpype.types import *
import os
import glob
import os
from datetime import datetime
import shutil

def get_latest_updated_directory(path):
    # 获取指定目录下的所有子目录
    subdirectories = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    
    # 如果没有子目录，则返回None
    if not subdirectories:
        return None

    # 获取最新更新的目录
    latest_subdir = max(subdirectories, key=os.path.getmtime)
    
    # 获取目录名
    latest_subdir_name = os.path.basename(latest_subdir)
    
    return latest_subdir_name

def run_dstap(network_name,javainit):
    """
    Runs the DSTAP partitioning algorithm with the specified arguments.

    Parameters:
        network_name (str): The name of the network to partition.
        sdda (str, optional): The SDDA value. Default is "2".
        partitioning (str, optional): The partitioning argument. Default is "partitioning".
    """
    # Define the paths to the directories containing the .class files and the .jar files
    class_path = r'C:\Users\DESH\Documents\Active\Box\Box Sync\2024Summer\1_DSTAP connect two system\Model File\DSTAP-master\build\classes'
    lib_path = r'C:\Users\DESH\Documents\Active\Box\Box Sync\2024Summer\1_DSTAP connect two system\Model File\DSTAP-master\build\classes\dstap\lib'

    # Create the classpath with all JAR files in the lib directory
    jars = glob.glob(os.path.join(lib_path, '*.jar'))
    classpath = [class_path] + jars

    # Start the JVM with the classpath
    if javainit:
        jpype.startJVM(classpath=classpath)

    # Import the Main class from the dstap.main package
    from dstap.main import Main

    # Define the arguments for the main method
    args_partitioning = ["partitioning", network_name, "sdda","2"]

    # Convert the args to Java array
    java_args = JArray(JString)(args_partitioning)

    # Call the main method of the Main class
    Main.main(java_args)
    folder_name=get_latest_updated_directory("C:\\Users\\DESH\\Documents\\Active\\Box\\Box Sync\\2024Summer\\1_DSTAP connect two system\\Model File\\DSTAP-master\\Networks\\SiouxFalls\\Inputs")
    # 定义源文件和目标文件路径
    source_file = os.path.join("C:\\Users\\DESH\\Documents\\Active\\Box\\Box Sync\\2024Summer\\1_DSTAP connect two system\\Model File\\DSTAP-master\\Networks\\SiouxFalls\\Inputs", "Parameters.txt")
    destination_file = os.path.join("C:\\Users\\DESH\\Documents\\Active\\Box\\Box Sync\\2024Summer\\1_DSTAP connect two system\\Model File\\DSTAP-master\\Networks\\SiouxFalls\\Inputs\\"+folder_name, "Parameters.txt")
    # 复制文件
    shutil.copy(source_file, destination_file)
    args_running = ["dstap_heuristic", "SiouxFalls", folder_name, "MEDIUM" , "1"]
    # Convert the args to Java array
    java_args = JArray(JString)(args_running)

    # Call the main method of the Main class
    Main.main(java_args)

    # Shut down the JVM
    #jpype.shutdownJVM()

#run_dstap("SiouxFalls")
