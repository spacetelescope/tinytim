/* File     :  system.c
 *
 * Contents :
 *      DeleteFile  :  Delete a file
 *      DefaultDir  :  Attach TINYTIM system directory
 *                       specifier to a filename.
 *
 * Author   :  John Krist (STScI)
 * Date     :  January 1992
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "tinytim.h"

/*---------------------------------------------------------------*/
void Delete_file( char *filename )
{
    unlink( filename );
}

/*----------------------------------------------------------------
* Default_dir :
*       Given a filename, attach the system-dependent directory
*       name for the TINYTIM system directory.
*
* Inputs :
*       string  :  Name of file in TINYTIM system directory.
*       filename : String to contain the fully specified filename.
*-----------------------------------------------------------------*/
void Default_dir( char *string, char *filename )
{

    char    *temp;
    if ( (temp = getenv("TINYTIM")) == NULL )
    {
        fprintf(stderr, "Warning: Environment variable TINYTIM undefined\n");
        if ( (temp = malloc(sizeof(char) * strlen(CONFIG_DEFAULT_DIR) + 1)) == NULL) {
            perror("Default_dir");
            exit(1);
        }
        fprintf(stderr, "Using default data directory: %s\n", CONFIG_DEFAULT_DIR);
        strcpy( temp, CONFIG_DEFAULT_DIR );
    }

    strcpy( filename, temp );
    strcat( filename, "/" );
    strcat( filename, string );

} /* Default_dir */
