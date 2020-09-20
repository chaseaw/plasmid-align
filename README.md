Prerequisites for using this pipeline
================================================================================
1. Linux-based command-line terminal
  *  Should work natively on Mac terminal. For Windows, it requires the linux subsystem.
  
* [WSL on Windows](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10)
  * Open Settings.
  
  *  Click on Apps.
    
  *  Scroll down to Related settings, and click on "Programs and Features".
    
  *  In the upper left corner, click on "Turn Windows features on or off".
    
  *  Make sure that "Windows Subsystem for Linux" is checked.
  
  *  Restart your computer.
  
  *  Download Ubuntu linux from the Windows Store (This is your terminal).
  
  *  First time using Ubuntu you will create a user name and password (Remember this!).

2. Git

*  From your terminal:

```
sudo apt update
sudo apt install git
```

3. IGV

Get IGV for your operating system (you don't need the command line version):
* [IGV Software Download](http://software.broadinstitute.org/software/igv/download)


Installation
================================================================================
From your home directory on the command line:

```
git clone https://github.com/chaseaw/plasmid-align.git
```

This will create a directory called plasmid-align that includes all necessary scripts.

Enter the new directory.

```
cd plasmid-align
```
Run the conda_setup script.
```
bash conda_setup.sh
```
You will be asked for your password. If you already have conda installed, this will only update package libraries and create a new environment.


Usage
================================================================================
