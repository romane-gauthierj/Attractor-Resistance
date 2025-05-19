import platform
import os
import stat
import shutil
import sys
import zipfile
from urllib.request import urlretrieve



def get_url(os_type, arch):
    base_url = "https://github.com/sysbio-curie/MaBoSS/releases/latest/download/"
    if os_type == 'linux':
        return base_url + "MaBoSS-linux64.zip"
    elif os_type == 'darwin':
        if arch == 'arm64':
            return base_url + "MaBoSS-osx-arm64.zip"
        elif arch == 'x86_64':
            return base_url + "MaBoSS-osx64.zip"
        else:
            print(f"Unsupported architecture: {os_type} : {arch}")
            sys.exit(1)
    elif os_type.startswith('win'):
        return base_url + "MaBoSS-win64.zip"
    else:
        print(f"Unsupported OS: {os_type}")
        sys.exit(1)
        
def ensure_directory_exists(path):
    """Create directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)
        
def download_and_extract(url, dest_path):
    """Download and extract a package from a URL."""
    # Ensure the destination directory exists
    ensure_directory_exists(dest_path)

    filename = url.split('/')[-1]
    file_path = os.path.join(dest_path, filename)

    # Download the file
    print(f"Downloading {filename} from {url}...")
    urlretrieve(url, file_path)

    # Extract the file
    print(f"Extracting {filename}...")
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(dest_path)

    for item in os.listdir(dest_path):
        st = os.stat(os.path.join(dest_path, item))
        os.chmod(os.path.join(dest_path, item), st.st_mode | stat.S_IEXEC)

    # Clean up
    os.remove(file_path)
    print(f"{filename} installed successfully.\n")
    
def is_installed(progname):
    return shutil.which(progname) is not None

def get_bin_path(os_type):
    
    if platform.system() == "Windows":
        bin_path = os.path.join(os.getenv("APPDATA"), "maboss", "bin")
    else:
        bin_path = os.path.join(os.path.expanduser("~"), ".local", "share", "maboss", "bin")

    if not os.path.exists(bin_path):
        os.makedirs(bin_path)
    
    if platform.system() == "Windows":
        os.environ["PATH"] = "%s;%s" % (bin_path, os.environ["PATH"])
    else:
        os.environ["PATH"] = "%s:%s" % (bin_path, os.environ["PATH"])
    return bin_path 

os_type = platform.system().lower()
arch = platform.machine().lower()

url = get_url(os_type, arch)
dest_path = get_bin_path(os_type)

if not is_installed("MaBoSS"):
    download_and_extract(url, dest_path)
    
print("MaBoSS is installed and ready to use.")
print(shutil.which("MaBoSS" if platform.system() != "Windows" else "MaBoSS.exe"))