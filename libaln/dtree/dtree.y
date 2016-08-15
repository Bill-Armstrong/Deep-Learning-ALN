%{

/* 
 * dtree.y
 * Dendronic Decisions Limited DTREE 1.1
 * syntax analysis
 */
                                      
%}

%union 
{                    
  double dbl;           /* numerical */
  long   n;             /* integer / flag */
  char*  psz;           /* pointer to null termintaed string */
}

/* reserved words */
%token  <n>     R_VERSION R_VARIABLES R_OUTPUT R_LINEARFORMS
%token  <n>     R_BLOCKS R_DTREE R_MIN R_MAX R_BLOCKREF

/* identifiers */
%token  <psz>   IDENT

/* special symbols */
%token  <n>     S_LESSEQUAL

/* constants */
%token  <dbl>   DBL_CONST
%token  <n>     INT_CONST

/* associativity */
%left '+' '-'


%start          DTREE_program    

%%

DTREE_program
  : version variables output linearforms blocks dtree
  ;

/*---  VERSION  ---*/

version
  : R_VERSION '=' DBL_CONST ';'
  ;

/*---  VARIABLES & OUTPUT  ---*/

variables
  : R_VARIABLES '=' dimension ';' varlist
  ;

dimension
  : INT_CONST
  ;

varlist
  : var 
  | varlist var 
  ;

var
  : varname ':' '[' var_minimum ',' var_maximum ']' ';'
  ;

varname
  : IDENT
  ;

var_minimum
  : DBL_CONST
  ;

var_maximum
  : DBL_CONST
  ;

output
  : R_OUTPUT '=' varname ';'
  ;

/*---  LINEARFORMS  ---*/

linearforms
  : R_LINEARFORMS '=' numlines ';' linelist
  ;

numlines
  : INT_CONST
  ;

linelist
  : line
  | linelist line
  ;

line
  : lineindex ':' lineexpr ';'
  ;

lineindex
  : INT_CONST
  ;

lineexpr
  : lineexpr '+' lineexpr
  | lineexpr '-' lineexpr
  | lineterm
  ;

lineterm
  : weight '(' centroidexpr ')'
  ;

centroidexpr
  : varname '+' centroid
  | varname '-' centroid
  ;

weight
  : DBL_CONST
  ;

centroid
  : DBL_CONST
  ;

/*---  BLOCKS  ---*/

blocks
  : R_BLOCKS '=' numblocks ';' blocklist
  ;

numblocks
  : INT_CONST
  ;

blocklist
  : block
  | blocklist block
  ;

block
  : blockindex ':' minmaxexpr ';'
  ;

blockindex
  : INT_CONST
  ;

minmaxexpr
  : lineindex
  | min_or_max '(' minmaxexprlist ')'
  ;

min_or_max
  : R_MIN
  | R_MAX
  ;

minmaxexprlist
  : minmaxexpr
  | minmaxexprlist ',' minmaxexpr
  ;

/*---  DTREE  ---*/

dtree
  : R_DTREE '=' numnodes ';' nodelist
  ;

numnodes
  : INT_CONST
  ;

nodelist
  : node
  | nodelist node
  ;

node
  : nodeindex ':' nodeexpr ';'
  ;

nodeindex
  : INT_CONST
  ;

nodeexpr
  : R_BLOCKREF blockindex
  | '(' varname S_LESSEQUAL threshold ')' '?' nodeindex ':' nodeindex
  ;

threshold
  : DBL_CONST
  ;

%%

