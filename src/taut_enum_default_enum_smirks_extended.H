//
// file taut_enum_default_enum_smirks_orig.H
// David Cosgrove
// AstraZeneca
// 2nd September 2011
//
// Contains the default SMIRKS definitions for the original tautomer
// enumeration as defined by Pete Kenny all those years ago, augmented by
// a far more aggressive set that do guanidinoguanidines and ketone-enols
// provided by (our) Andrew Leach.

#ifndef TAUT_ENUM_DEFAULT_ENUM_EXTENDED_H
#define TAUT_ENUM_DEFAULT_ENUM_EXTENDED_H

#include <string>

namespace DACLIB {

const std::string ENUM_SMIRKS_EXTENDED = "\n\
# May 2012\n\
# The set described below was getting very complex and causing more problems than it solved.\n\
# The set has therefore been cut back to rules that should be more straightforward.\n\
# There are a set of exceptions for some aromatic rings that are specific problem cases for the \n\
# various aromaticity models.\n\
# Aromatic 7 membered rings are also a problem that I fear I have not properly fixed.\n\
#\n\
# November 2011\n\
# This is a set of SMIRKS from two sources and intended to achieve complete enumeration of tautomers\n\
# Previously, the protocol for achieving this applied a set of rules for standardising, some of which\n\
# actually involved selecting an arbitrary tautomer for certain tautomer types (eg phenol vs aromatic\n\
# amide).  Standardising was followed by a limited enumeration.  The SMARTS for this limited enumeration \n\
# are at the top of this file and are called ENUM_#.  At the bottom, the list of standardising rules that involve tautomer\n\
# changes are included and have names ending _ENUMFWD.  Paired with most of these are rules for moving \n\
# the protons back to generate the other tautomer and these rules end _ENUMBWD.  Some of these\n\
# required some new vector binding definitions which have been added in the accompanying vb file.\n\
# Some of the rules such as N=C-OH in aliphatic systems probably do not want to be reversed and\n\
# so do not have a pair and are not in this set.  \n\
#\n\
# ALIPHATIC SYSTEMS\n\
#\n\
[N;X2:1]=[C:2]-[O,S;H;X2:3][H:4]>>[N:1]([H:4])-[C:2]=[*:3] ENUM_ALI_1\n\
# This one has been expanded to take in imine vs enamine \n\
# Fixed to not make allene like structures\n\
[H:4][O,S,N;!H0:1]-[C,N:2]([#1,C,c:5])=[C:3]>>[*:1]=[*:2]([*:5])-[*:3][H:4] ENUM_ALI_2F\n\
[O,S,N!$(N(=*)(=*)*)!+:1]=[C,N:2]([#1,C,c:5])-[C;!H0!$(C=*):3][H:4]>>[H:4][*:1]-[*:2]([*:5])=[*:3] ENUM_ALI_2B\n\
[H:4][N;!H0:1]-[C,N:2]=[C,N!$(N(=*)(=*)*)!+:3]>>[*:1]=[*:2]-[*:3][H:4] ENUM_ALI_3F\n\
[N!$(N(=*)(=*)*)!+:1]=[C,N:2]-[C,N;!H0!$([C,N]=*):3][H:4]>>[H:4][*:1]-[*:2]=[*:3] ENUM_ALI_3B\n\
#\n\
[O,S&X2,N&X3;!+:1]1[NX2:2]=[C:3]-[C;H,H2:4]([H:7])[C:5]1=[O,S:6]>>[O,S&X2,N&X3;!+:1]1[N:2]([H:7])-[C:3]=[C:4][C:5]1=[O,S:6] ENUM_ALI_4F\n\
[O,S&X2,N&X3;!+:1]1[N!H0:2]([H:7])-[C:3]=[C:4][C:5]1=[O,S:6]>>[O,S&X2,N&X3;!+:1]1[NX2:2]=[C:3]-[C;H,H2:4]([H:7])[C:5]1=[O,S:6] ENUM_ALI_4B\n\
#\n\
[NH:1]=[C,N:2]-[NH:3]([H:5])-[*:4]>>[H:5][NH:1]-[*:2]=[NH:3]-[*:4] ENUM_ALI_5F\n\
[H:5][NH:1]-[C,N:2]=[NH:3]-[*:4]>>[NH:1]=[*:2]-[NH:3]([H:5])-[*:4] ENUM_ALI_5B\n\
#AGL Extras\n\
[OH:1]([H:4])-[N:2]=[A!$(N(=*)(=*)*):3]>>[O:1]=[N:2]-[*:3][H:4] ENUM_ALI_6F\n\
[O:1]=[N:2]-[A!$(N(=*)(=*)*):3][H:4]>>[OH:1]([H:4])-[N:2]=[*:3] ENUM_ALI_6B\n\
#\n\
# AROMATIC SYSTEMS\n\
#\n\
#nH to n\n\
[H:4][n;H1;X3;!+:1]:[$AR:2]:[$NPY:3]>>[n:1]:[*:2]:[nH:3][H:4] ENUM_AROM_7_1\n\
[H:6][n;H1;X3;!+:1]:[$AR:2]:[$AR:3]:[$AR:4]:[$NPYNoImid:5]>>[n:1]:[*:2]:[*:3]:[*:4]:[nH:5][H:6] ENUM_AROM_7_3\n\
[H:8][n;H1;X3;!+:1]:[$AR:2]:[$AR:3]:[$AR:4]:[$AR:5]:[$AR:6]:[$NPYNoImid:7]>>[n:1]:[*:2]:[*:3]:[*:4]:[*:5]:[*:6]:[nH:7][H:8] ENUM_AROM_7_5\n\
[H:4][$OS:1]-[c:2]:[$NPY:3]>>[*:1]=[*:2]:[*:3][H:4] ENUM_AROM_8_1\n\
[H:6][$OS:1]-[c:2]:[$AR:3]:[$AR:4]:[$NPY:5]>>[*:1]=[*:2]:[*:3]:[*:4]:[*:5][H:6] ENUM_AROM_8_3\n\
[H:8][$OS:1]-[c:2]:[$AR:3]:[$AR:4]:[$AR:5]:[$AR:6]:[$NPY:7]>>[*:1]=[*:2]:[*:3]:[*:4]:[*:5]:[*:6]:[*:7][H:8] ENUM_AROM_8_5\n\
# NH to O or S \n\
# n changed to #7 to pick up one some problem cases and : bond definition removed\n\
[H:4][#7;H1;X3;!+:1]:[c:2]=[$REV_OS:3]>>[*:1]:[c:2]-[*:3][H:4] ENUM_AROM_9_0\n\
[H:6][#7;H1;X3;!+:1]:[$AR:2]:[$AR:3]:[c:4]=[$REV_OS:5]>>[*:1]:[*:2]:[*:3]:[c:4]-[*:5][H:6] ENUM_AROM_9_2\n\
[H:8][#7;H1;X3;!+:1]:[$AR:2]:[$AR:3]:[$AR:4]:[$AR:5]:[c:6]=[$REV_OS:7]>>[*:1]:[*:2]:[*:3]:[*:4]:[*:5]:[c:6]-[*:7][H:8] ENUM_AROM_9_4\n\
#OH or SH to O or S\n\
[H:6][$OS,$EXO_NH:1]-[c:2]:[a:3]:[c:4]=[$OS,$EXO_NDB:5]>>[*:1]=[c:2]:[*:3]:[c:4]-[*:5][H:6] ENUM_AROM_10_1\n\
[H:8][$OS,$EXO_NH:1]-[c:2]:[a:3]:[a:4]:[a:5]:[c:6]=[$OS,$EXO_NDB:7]>>[*:1]=[c:2]:[*:3]:[*:4]:[*:5]:[c:6]-[*:7][H:8] ENUM_AROM_10_3\n\
[H:10][$OS,$EXO_NH:1]-[c:2]:[a:3]:[a:4]:[a:5]:[a:6]:[a:7]:[c:8]=[$OS,$EXO_NDB:9]>>[*:1]=[c:2]:[*:3]:[*:4]:[*:5]:[*:6]:[*:7]:[c:8]-[*:9][H:10] ENUM_AROM_10_5\n\
#dealing with weird 7 membered aromatics\n\
[n;H;$(n1aaaaaa1):1]([H:4]):[c:2]=[N:3]>>[*:1]:[*:2]-[*:3][H:4] STAND_84_ENUMFWD\n\
[N;$([#7]1******1):1]=[C;$([#6]1******1):2]-[N;!H0:3][H:4]>>[H:4][*:1]-[*:2]=[*:3] STAND_85_ENUMFWD\n\
[n;$([#7]1******1):1]:[c;$([#6]1******1):2]-[N;!H0:3][H:4]>>[H:4][*:1]:[*:2]=[*:3] STAND_86_ENUMFWD\n\
# And some extended systems that the aromatic model messes up\n\
[nH:1]1([H:35]):[c:2]2[c:3]3:[c:4](=[O:5])[c:6]4[c:7][c:8][c:9][c:10][c:11]4[c:12](=[O:13]):[c:14]3[c:15][c:16][c:17]2:[nH:18]([H:36])[c:19]5[c:20]6[c:21](=[O:22])[c:23]7[c:24][c:25][c:26][c:27][c:28]7[c:29](=[O:30])[c:31]6[c:32][c:33][c:34]51>>[N:1]1=[c:2]2[c:3]3=[C:4](-[O:5][H:35])[c:6]4[c:7][c:8][c:9][c:10][c:11]4[C:12](-[O:13][H:36])=[c:14]3[c:15][c:16][c:17]2=[N:18][c:19]5[c:20]6[c:21](=[O:22])[c:23]7[c:24][c:25][c:26][c:27][c:28]7[c:29](=[O:30])[c:31]6[c:32][c:33][c:34]51 ENUM_AROM_AZ1F\n\
[N:1]1=[c:2]2[c:3]3=[C:4](-[O:5][H:35])[c:6]4[c:7][c:8][c:9][c:10][c:11]4[C:12](-[O:13][H:36])=[c:14]3[c:15][c:16][c:17]2=[N:18][c:19]5[c:20]6[c:21](=[O:22])[c:23]7[c:24][c:25][c:26][c:27][c:28]7[c:29](=[O:30])[c:31]6[c:32][c:33][c:34]51>>[nH:1]1([H:35]):[c:2]2[c:3]3:[c:4](=[O:5])[c:6]4[c:7][c:8][c:9][c:10][c:11]4[c:12](=[O:13]):[c:14]3[c:15][c:16][c:17]2:[nH:18]([H:36])[c:19]5[c:20]6[c:21](=[O:22])[c:23]7[c:24][c:25][c:26][c:27][c:28]7[c:29](=[O:30])[c:31]6[c:32][c:33][c:34]51 ENUM_AROM_AZ1B\n";

}

#endif // TAUT_ENUM_DEFAULT_ENUM_EXTENDED_H
