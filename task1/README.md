# Welcome to BIOL 469, Genomics

This github repository contains a series of tasks designed to teach students:

  * Command-line linux
  * Genomics & bioinformatics data analysis


## Getting Started

Before you can begin with the coding exercises, you must have access to a linux machine.
You can either use your own local system or a remote VM that has been set up for you.

### Accessing the remote VM

If you have been given access to a remote VM, you can access it either through the url that has been given to you, or directly through the Terminal like this:

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

## Once you are logged on | Learning Unix

If you have logged in correctly, you should see a welcome screen.

You are now in a unix command-line environment. For more information and instructions:

* [CodeAcademy](https://www.codecademy.com/learn/learn-the-command-line) - Learn the command line
* [TeachingUnix](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Another Unix Tutorial
* [BasicCommands](http://mally.stanford.edu/~sr/computing/basic-unix.html) - List of common commands

## A mini-primer

To make the folder "task1", type

```
mkdir task1
```

Change directory into 'task1' folder

```
cd task1
```

And now create a file called file.txt

```
>file.txt
```

Delete the file

```
rm file.txt
```

Move back to the previous folder

```
cd ..
```

And delete the folder 'task1'

```
rmdir task1
```


Congratulations, you're done. You are now ready for task2.


