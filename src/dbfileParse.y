%{

#include <stdio.h>
#include "dbimage.h"

#undef yywrap
int     yywrap          ( void );
int     yyparse(void);

void where_am_i();
int yylex();
FILE *yyin;

%}

%union {
  char *string;
}

%token IMAGEPATH
%token  <string> STRING

%type <string> image_path
%type <string> image_pathes

%%

config : definitions
;

definitions : definition
            | definitions definition

definition : IMAGEPATH '{' image_pathes '}'
;

image_pathes : image_path
        | image_pathes  image_path
;


image_path : STRING ':' STRING { DbConfigAddImagePath($3,$1); $$ = $1;}
         | image_path ',' STRING {  DbConfigAddImagePath($3,$1); $$ = $1;}
;
%%

void yyerror(char *Message)
{
printf("%s\n",Message);
where_am_i();
}

int yywrap( void ) {return( 1 );}

int DbConfigFileParse(const char * FileName)
{
yyin  = fopen(FileName,"r");
if (!yyin) 
  {
  printf(" DbConfigFileParse  : cannot open %s\n", FileName);
  return 0;
  }
return (yyparse() == 0);
}

