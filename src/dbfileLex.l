%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include  "dbfileParse.h"

static char *StringNew(char *Original)
{
char *toto = (char*) malloc(strlen(Original)+2);
strcpy(toto,Original);
return toto;
}

int yywrap ();
int yylook ();
int yyback ();


%}

%option yylineno


sep [\{\}:,]
non_sep [^\{\}:,\n\t ]
space [ \t\n]
non_space [^  \t\n]
%%
^#.*$              {;}
"CatalogPath"{space}* {return CATALOGPATH;}
"ImagePath"{space}* {return IMAGEPATH;}
\. {yylval.string = StringNew(yytext); return STRING;}
{non_space}+ {if (strlen(yytext) == 1)
                {
                   return yytext[0]; 
                }
              else 
                 {
                   yylval.string = StringNew(yytext); 
                   return STRING;
                 }
              }
{space}+   {;}
\n			   {;}
%%

/* last rule may be '.|\n' */

void where_am_i() 
/* to tell the (poor) user where the scanner is when a syntax error occurs */
{
printf(" line number %d, reading %s \n",yylineno, yytext);
}

/* int yywrap() { return 1;} */
