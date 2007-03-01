#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXVAR 400

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
#define MAX_LENGTH 30
char **tags = (char **) calloc(MAXVAR,sizeof(char*));
int i, dim = 0;
char line[512];

while (fgets(line,512,file))
  {
  if (line[0] == '@') continue;
  char w1[512], w2[512];
  if (line[0] != '#')
    {
    printf(" header should end by '# end'\n"); return NULL;
    }
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

char *split_line(char *Line, float *X, const int Dim, int *nread)
{
  const char *p1;
  char* p2;
  int i;
  *nread = 0;
  if (strlen(Line) <= 1) return Line;
  memset(X,0,Dim*sizeof(X[0]));
  p1 = Line;
  for (i=0; i< Dim; ++i)
    {
    float value = strtod(p1,&p2);
    if (p2 == p1) break;
    X[i] = value;
    p1 = p2;
    (*nread)++;
    }
  return p2;
}



int read_data(FILE *File, int Dim, void (Processor)(int*,float*), int *ProcData)
{
float x[MAXVAR];
char line[8192];
 int count = 0;
int miss = 0; int more = 0;

 if (Dim>MAXVAR)
   {
     Dim = MAXVAR;
   }
while (fgets(line,8192,File))
  {
  char *left_over;
  int nread;
  if (strlen(line) <= 1) continue;
  if (line[0] == '#' || line[0] == '@') continue;
  left_over = split_line(line,x,Dim,&nread);
  if (nread < Dim) miss++;
  if (atof(left_over)) more++;
  if (Processor) Processor(ProcData,x);
  count ++;
  }
if (miss || more)
  {
  printf(" when reading colums, we missed items on %d rows and had more on %d rows \n",miss, more);
  }
 return count;
}

