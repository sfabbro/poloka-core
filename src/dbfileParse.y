%{

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dbimage.h"

#undef yywrap
int     yywrap          ( void );
int     yyparse(void);

void where_am_i();
int yylex();
FILE *yyin;

static char* StringCopy(const char *In)
{
  char *toto = (char *) malloc(strlen(In)+2);
  strcpy(toto,In);
  return toto;
}

static char* StringCat(const char *In1, const char *In2)
{
  char *toto = (char *) malloc(strlen(In1)+strlen(In2)+3);
  sprintf(toto,"%s %s",In1,In2);
  return toto;
}


void yyerror(char *Message)
{
printf("%s\n",Message);
where_am_i();
}


%}

%union {
  char *string;
}

%token IMAGEPATH
%token  <string> STRING

%type <string> image_path
%type <string> image_pathes

%token CATALOGPATH
%type <string> catalog_path
%type <string> catalog_pathes

%token IMAGENAMES
%type <string> names
%%

config : definitions
;

definitions : definition
            | definitions definition
;

definition : IMAGEPATH '{' image_pathes '}'
            | CATALOGPATH '{' catalog_pathes '}'
            | IMAGENAMES '{' image_names '}'
;

image_pathes : image_path
        | image_pathes  image_path
;

catalog_pathes : catalog_path
        | catalog_pathes catalog_path
;


catalog_path : STRING {DbConfigAddCatalogPath($1); $$ = $1;}
;

image_path : STRING ':' STRING { DbConfigAddImagePath($3,$1); $$ = $1;}
         | image_path ',' STRING {  DbConfigAddImagePath($3,$1); $$ = $1;}
;


image_names : image_name_block
         | image_names image_name_block
;

image_name_block : STRING '{' names '}' 
                {DbConfigAddNewImageNames($1,$3); free($3);}
            | '{' names '}' {DbConfigAddNewImageNames("",$2); free($2);}
;

names : STRING { $$ = StringCopy($1);}
        | names STRING { $$ = StringCat($1,$2); free($1);}
;
%%


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

