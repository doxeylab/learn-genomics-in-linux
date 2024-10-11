# Installing Artemis/Java

Artemis Installation Instructions

## Part 1 - Install the Java Development Kit (JDK):
1.	Go to the [Temurin by Adoptium website](https://adoptium.net/en-GB/temurin/releases/?version=11&os=windows&arch=x64&package=jdk) and select your operating system version (most systems run x64 bit).
2.	Select the "**JDK 11-LTS**" version from the dropdown menu. 
-	For Windows: Download the **.msi** file and run it once downloaded.
-	For MacOS: Download the **.pkg** file and run it once downloaded.

## Part 2 - Install Artemis Tools:
1.	Go to the [Sanger Pathogens website](https://sanger-pathogens.github.io/Artemis/) and scroll down the page until you reach the **"Software Availability"** section.
2.	Under the **"Download"** heading, select your operating system version from the list and download the file. 
-	For Windows: Download the **.zip** file.
-	For MacOS: Download the **.dmg** file.
3.	Unzip the file to an appropriate local directory and an Artemis folder will be created containing the tools. 
-	For Windows: Recommended to unzip the file to the **"C:\"** drive.
-	For MacOS: Recommended to unzip the file to the "**Applications**" folder.
The first time you run one of the tools, it will ask to "Set Working Directory". Ask your prof what to set this to as they probably have some files to get data from – the working directory should be the one that contains those.

#### Helpful resource:
[The Artemis Manual](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.html)

## Additional Steps For MacOS:
You may run into the error message "This application requires that Java 9 or later be installed on your computer. Please download and install the latest version of Java from www.java.com and try again." when opening any of the tools. If so, please follow the instructions below to open Artemis.
1.	Open the "**Artemis**" folder.
2.	Right click on the "**Artemis**" icon and select "Show Package Content".
3.	Double-click on the "**Contents**" folder.
4.	Double-click on the "art" executable file (it has the black terminal icon).
5.	A new window will open and you will be prompted to set the working directory.
6.	Once confirmed it is working, you can follow the additional steps below to create a shortcut for each of the Artemis tools. 
-	Go through steps 1-3 again until the "**Contents**" folder is open.
-	Right-click on the "art" executable file and from the dropdown menu, click "**Make Alias**".
-	A new shortcut called "**art alias**" will be made in the folder. You can drag this file to your Desktop for example for easy access.
-	From the Desktop (or wherever you put the file), double-click on "**art alias**" and the program should open without errors.
-	You can repeat these steps for ACT, BamView, and Circular-Plot.

