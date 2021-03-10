Analysis of Nextera-NGS sequenced plasmids
================================================================================
Description: This directory will enable easy analysis of plasmids sequenced by Nextera XT library prep.

There are a number of ways to do this but this one should be comprehensive and relatively painless.

[Go To Installation](#installation)

[Go To Usage](#usage)

Prerequisites for using this pipeline
================================================================================
**1.** **Linux-based command-line terminal**
  *  Should work natively on Mac terminal. For Windows, it requires the linux subsystem.
  
* [WSL on Windows](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10)
  * Open Settings.
  
  *  Click on Apps.
    
  *  Under Related settings, click on "Programs and Features".
    
  *  In the upper left corner, click on "Turn Windows features on or off".
    
  *  Make sure that "Windows Subsystem for Linux" is checked.
  
  *  Restart your computer.
  
  *  Download Ubuntu from the Microsoft Store App (Ubuntu will be your terminal).
  
  *  First time using Ubuntu you will create a user name and password (Remember this!).

**2.** **Git**

*  From your terminal:

```
sudo apt update
sudo apt install git
```

**3.** **IGV**

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
You may be asked for your password. If conda is already installed, this will only update package libraries and create a new environment.

Sometimes conda does not properly initialize. If an error appears, usually running ```bash conda_setup.sh``` again will fix it. 

To make sure everything is working properly, run the plasmid-align.sh script in place on included test data:
```
source plasmid-align.sh Test_Data*
```
You should see the script running. Test_Data directories should now both have populated de_novo and ref_alignment folders.

The Test_Data_1 will output a perfect alignment to reference, while Test_Data_2 will show examples of mutations, deletions, and insertions by IGV.

Usage
================================================================================
make a new directory to put reads into, and change into your new directory.

```
mkdir <Today's Date>_Plasmid_Sequencing
cd <Today's Date>_Plasmid_Sequencing
```

retrieve your files from the link provided by the SIC.

```
wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://htcf.wustl.edu/files/<Your Directory>/
```
The --cut-dirs=2 option refers to the number of directories between the .edu site and your files (which is 2 in the example).

copy the Sample-folders.sh and plasmid-align.sh files to this new directory containing plasmid sequencing paired-end reads still in Illumina format fastq.gz.

An easy way to copy or move files between these folders in windows is to open the linux subsystem in explorer

```
explorer.exe .
```

Alternatively, you can use the copy command in bash to move the run files to your folder.

From Your Directory:
```
cp ../Sample-folders.sh ./Sample-folders.sh
cp ../plasmid-align.sh ./plasmid-align.sh
```
This assumes that you made your directory in the plasmid-align folder (thus ../ or up one directory) and that you are copying to your current directory (./)

From Your Directory, run the Sample-folders.sh script to unzip and separate R1 and R2 samples into their own folder.

```
bash Sample-folders.sh
```
Each Sample folder needs a reference fasta sequence file (.fa) if you want to align to a known sequence.

Here I find that using the ```explorer.exe .``` command is very helpful for dropping fasta files into the Sample Folders.

After adding each fasta file, from the directory with Sample folders, run plasmid-align.sh:

```
source plasmid-align.sh Sample_*
```

The * wildcard ensures all Sample folders will be analyzed. Could be run on a single folder as well.

output files are final contigs .fa in the de_novo folder and the bam files in the ref_alignment folder.

The longest contigs in the final contigs .fa file are usually the entire plasmid sequence, but not always perfect with repetitive regions.

Run IGV loading the reference .fa as genome and your new .bam as a track to identify any sequence changes.
