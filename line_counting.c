#include<stdio.h>
#include<stdlib.h>

int read_file_line(char a[])
{
    FILE *fp;
    char ch;
    fp = fopen(a,"r");
    int lines=1;

    if (fp == NULL)
    {
        printf("\nError while opening the file\n");
        return(0);
    }

    // EOF = End of file
    while((ch=fgetc(fp))!=EOF)
    {
        if(ch == '\n')
        {
            lines++;
        }
    }


    return lines;
}