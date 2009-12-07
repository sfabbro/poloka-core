#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXVAR 4096

static char *StringNew(char *Source) /* prototype */
{
char *toto;
toto = (char *) calloc(1,strlen(Source)+1);
strcpy(toto,Source);
return toto;
}

static char * CutExtension(const char *Name)
{
static char result[512];
 char *dot;
strcpy(result,Name);
dot = strrchr(result, '.');
if (dot) *dot = '\0';
return result;
}

char** decode_tags(FILE *file, int *Dim)
{
#define MAX_LINE_LENGTH 16384
char **tags = (char **) calloc(MAXVAR,sizeof(char*));
int i, dim = 0;
char line[MAX_LINE_LENGTH];

while (fgets(line,MAX_LINE_LENGTH,file))
  {
  if (line[0] == '@') continue;
  if (line[0] != '#')
    {
    printf(" header should end by '# end'\n"); return NULL;
    }
  char w1[MAX_LINE_LENGTH], w2[MAX_LINE_LENGTH];
  if (sscanf(line+1,"%s",w1) == 1 && strcmp(w1,"end") == 0) break;
  if (dim == MAXVAR) 
    { printf(" tuple truncated to %d variables\n",dim); continue;} 
  if (sscanf(line+1," %[^: ] %s",w1,w2)==2  && w2[0] == ':' )
    {
      tags[dim] = StringNew(w1);(dim)++;
      continue;
    }
  } /* end of input read */
*Dim = dim;
if (dim == 0) return 0;
for (i =0; i<dim; ++i) printf(" %d %s\n",i,tags[i]);
return tags;
}

char *split_line(char *Line, float *X, const int Dim, int *nread, int *nbad)
{
  const char *p1;
  char dummy[MAX_LINE_LENGTH];
  char* p2;
  int i;
  *nread = 0;
  *nbad =0;
  if (strlen(Line) <= 1) return Line;
  memset(X,0,Dim*sizeof(X[0]));
  p1 = Line;
  for (i=0; i< Dim; ++i)
    {
    float value = strtod(p1,&p2);
    if (p2 == p1) /* try to read a bunch of chars to go on */
      {
	int nread=0;
	sscanf(p1,"%s%n",dummy,&nread);
	if (nread == 0) break; // means end of line
	(*nbad)++;
	p1 += nread;
	X[i] = 1e30;
      }
    else
      {
	X[i] = value;
	p1 = p2;
      }
    (*nread)++;
    }
  return p2;
}


int read_data(FILE *File, int Dim, void (Processor)(int*,float*), 
	      int *ProcData, bool print_bad_lines)
{
  char line[MAX_LINE_LENGTH];
  int count = 0;
  int miss = 0; int more = 0;
  int bad = 0;
  float x[MAXVAR];

  if (Dim>MAXVAR)
    {
      Dim = MAXVAR;
    }
  while (fgets(line,MAX_LINE_LENGTH,File))
    {
      char *left_over;
      int nread;
      int nbad;
      if (line[0] == '#' || line[0] == '@') continue;
      size_t len=strlen(line);
      if (len <= 1) continue;
      if (line[len-1] == '\n') line[len-1] = '\0'; 
      left_over = split_line(line,x,Dim,&nread, &nbad);
      bool error = false;
      if (nread < Dim) {miss++; error=true;}
      if (atof(left_over)) {more++; error = true;}
      if (nbad) {bad++; error = true;}
      if (print_bad_lines && error)
	printf(" bad line :\n%s\n",line);
      if (Processor) Processor(ProcData,x);
      count ++;
    }
  if (miss )
    printf(" when reading colums, we missed items on %d rows \n",miss);
  if (more)
    printf(" when reading colums, we had too many items on %d rows \n",more);
  if (bad)
    printf(" when reading colums, we had bad conversions on %d rows\n",bad);
  return count;
}


