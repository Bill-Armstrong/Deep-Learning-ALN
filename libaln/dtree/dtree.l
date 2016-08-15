%{

/* 
 * dtree.l
 * Dendronic Decisions Limited DTREE 1.1 
 * lexer
 */
                                        
#include "ytab.h"

%}
            
comment   "//".*\n        
ws        [ \t\r\n]+
digit     [0-9]
letter    [a-zA-Z]
integer   -?{digit}+
float     -?({integer}\.{integer})|({integer}(\.{integer})?[eE][+-]?{integer})
ident     {letter}({letter}|_|{digit})*
   
%%

{ws}          ;
{comment}     ;                      
{float}       { yylval.dbl = atof(yytext); return DBL_CONST; }
{integer}     { yylval.n = atol(yytext); return INT_CONST; }
"VERSION"     { yylval.n = R_VERSION; return R_VERSION; }
"VARIABLES"   { yylval.n = R_VARIABLES; return R_VARIABLES; }
"OUTPUT"      { yylval.n = R_OUTPUT; return R_OUTPUT; }
"LINEARFORMS" { yylval.n = R_LINEARFORMS; return R_LINEARFORMS; }
"BLOCKS"      { yylval.n = R_BLOCKS; return R_BLOCKS; }
"block"       { yylval.n = R_BLOCKREF; return R_BLOCKREF; }
"DTREE"       { yylval.n = R_DTREE; return R_DTREE; }
"MIN"         { yylval.n = R_MIN; return R_MIN; }
"MAX"         { yylval.n = R_MAX; return R_MAX; }
"<="          { yylval.n = S_LESSEQUAL; return S_LESSEQUAL; }
{ident}       { yylval.psz = yytext; return IDENT; /* make sure to *copy* this string if needed! */ } 
.             { yylval.n = *yytext; return *yytext; }

%%
