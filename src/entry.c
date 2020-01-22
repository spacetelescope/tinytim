/*  File : entry.c
 *  Contents :
 *     Init_entry_list   : Initialize counter for parameter buffer.
 *     Delete_entry_list : Deallocate the parameter buffer.
 *     Write_entry_list  : Write current parameter buffer to a file.
 *     Store_entry      : Create an entry in the parameter buffer.
 *     Get_string       : Read a line from a file.  Does NOT store line.
 *     Get_entry        : Read a line from a file and store it in the buffer.
 *
 *  Written by John Krist, June 1993.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "tinytim.h"

#ifndef MACOSX
#include <malloc.h>
#endif
#include <memory.h>

static  char  *Entry_list[1000];
static  int   Entry_number;
static   char  Entry[MAX_STRING];

/*---------------------------------------------------------------------------*/
void Init_entry_list( void )
{
    Entry_number = 0;
}

/*---------------------------------------------------------------------------*/
void Delete_entry_list( void )
{
    int     i;

    for ( i = 0; i < Entry_number; ++i )
        free( Entry_list[i] );
}

/*---------------------------------------------------------------------------*/
void Write_entry_list( FILE *file )
{
    int     i;

    for ( i = 0; i < Entry_number; ++i )
        fputs( Entry_list[i], file );
}

/*---------------------------------------------------------------------------*/
void Store_entry( char *entry )
{
    Entry_list[Entry_number] = (char *)malloc( MAX_STRING );
    strcpy( Entry_list[Entry_number], entry );
    ++Entry_number;
}

/*--------------------------------------------------------------------------
* Routine : Get_string
* Purpose : Read a line from a table, ignoring comment (#) lines.  Unlike
*           Get_entry(), it does not store the line read into a list.
*--------------------------------------------------------------------------*/
char *Get_string( FILE *file )
{
    char    *n;

    if ( fgets( Entry, MAX_STRING - 1, file ) == NULL )
        return( NULL );

    while ( Entry[0] == '#' )
        if ( fgets( Entry, MAX_STRING - 1, file ) == NULL )
            return( NULL );

    /* Remove any comments at the end of the line by placing a \0  *
     * after the last parameter value on the line.                 */

    n = strchr( Entry, '#' );
    if ( n != NULL )
    {
        --n;

        /* Work backwards to the last visible character */

        while ( !isgraph(*n) )
            *n-- = '\0';
    }

    return( Entry );

} /* Get_string */

/*--------------------------------------------------------------------------
* Routine : Get_entry
* Purpose : Read a line from a table, ignoring comment (#) lines. Stores
*           all lines read into a list.
*--------------------------------------------------------------------------*/
char *Get_entry( FILE *file )
{
    char    *n;

    if ( fgets( Entry, MAX_STRING - 1, file ) == NULL )
    {
        fprintf(stderr, "Get_entry : Unexpected EOF (line=%d)\n",
            Entry_number );
        exit(2);
    }

    Store_entry( Entry );

    while ( Entry[0] == '#' )
    {
        if ( fgets( Entry, MAX_STRING - 1, file ) == NULL )
        {
            fprintf(stderr, "Get_entry : Unexpected EOF (line=%d)\n",
                Entry_number );
            exit(2);
        }
        Store_entry( Entry );
    }

    /* Remove any comments at the end of the line by placing a \0  *
     * after the last parameter value on the line.                 */

    n = strchr( Entry, '#' );
    if ( n != NULL )
    {
        --n;

        /* Work backwards to the last visible character */

        while ( !isgraph(*n) )
            *n-- = '\0';
    }

    return( Entry );

} /* Get_entry */

/* A gross hack to implement the focus adjustment approach documented
 * in the user's guide as expediently as possible.  Whatever z4 is,
 * add waves of secondary mirror despace to it.
 */
void Hack_z4_entry(float waves_secondary_mirror_despace)
{
    int i;
    for(i=0; i<Entry_number; i++) {
        if (strstr(Entry_list[i], "# Z4") || strstr(Entry_list[i], "# z4")) {
            float z4;
            char comment[MAX_STRING];
            sscanf(Entry_list[i], "%f # %[^\n]", &z4, comment);
            z4 += waves_secondary_mirror_despace;
            sprintf(Entry_list[i], "%f # %s\n", z4, comment);
        }
    }
}
