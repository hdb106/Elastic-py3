#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

int stoi(char *s, int *iar)
{
  int j=0,k,l=0,n=0;
  char str[128];

  while (1)
  {
    while ( s[j] == ' ' || s[j] == '\t' || s[j]==',') j++;
    if( s[j] == '\0' || s[j] == '/' || s[j] == '\n' ) return n;
    k=0;
    while ( isdigit(s[j]) || s[j]=='+' || s[j]=='-')  str[k++]=s[j++];
    str[k]='\0';
    iar[l++]=atoi(str);
    n++;
  }
}

int stod(char *s, double *dar)
{
  int j=0,k,l=0,n=0;
  char str[128];

  while (1)
  {
    while ( s[j] == ' ' || s[j] == '\t' || s[j]==',') j++;
    if( s[j] == '\0' || s[j] == '/' || s[j] == '\n' ) return n;
    k=0;
    while ( isdigit(s[j]) || s[j]=='.' || tolower(s[j])=='e' ||
            s[j]=='+' || s[j]=='-')  str[k++]=s[j++];
    str[k]='\0';
    dar[l++]=atof(str);
    n++;
  }
}

int stos(char *s, char *s_out)
{
  int j=0,k;

  while ( s[j] == ' ' || s[j] == '\t' || s[j]==',') j++;
  if( s[j] == '\0' || s[j] == '/' || s[j] == '\n' ) return 0;
  k=0;
  while ( isdigit(s[j]) || isalpha(s[j]))  s_out[k++]=s[j++];
  s_out[k]='\0';
  return 1;
}

void print_dbl(char *text, int n, double *d)
{
   int i;
   if( text !=NULL ) printf("%s",text);
   for(i=0; i<n; i++)
     printf(" % .8f",d[i]);
   printf("\n");
   return;
}
