# Task1 - accessing and navigating the Linux command-line

This task will introduce students to the Linux command-line.

## Getting Started

Before you can begin with the coding exercises, you must have access to a linux machine.
You can either use your own local system or a remote VM that has been set up for you.

### Accessing the remote VM

If you have been given access to a remote VM, you can access it either through the url that has been given to you, or directly through the Terminal like this:

```
ssh -i /path/to/your/.ssh/publickey yourUserName@remoteIP
```

When you are done, you can leave your session by typing...

```
exit
```

## Once you are logged on | Learning Unix

If you have logged in correctly, you should see a welcome screen.

You are now in a unix command-line environment. For more information and instructions:

* [CodeAcademy](https://www.codecademy.com/learn/learn-the-command-line) - Learn the command line
* [TeachingUnix](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Another Unix Tutorial
* [BasicCommands](http://mally.stanford.edu/~sr/computing/basic-unix.html) - List of common commands

## A unix primer

### Navigating files/folders

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

Open the file with `nano`

```
nano file.txt
```

Now enter a few lines of text, type 'ctrl-o' and then 'ctrl-x' to save and exit

To print the contents of your file type

```
cat file.txt
```

To count the number of words and lines in your file

```
wc file.txt #words
wc -l file.txt #lines
```

To copy the file to a new file

```cp file.txt newfile.txt```

To combine both files together into a third file

```
cat file.txt newfile.txt > thirdfile.txt
```

The '>' redirects the output of the commands on the left of it to a file specified on the right.

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


