
from developer import DevBase


#
# pre-made libraries to leverage
#   spglib 
#     http://spglib.sourceforge.net/
#     -Spglib is a C library for finding and handling crystal symmetries.
#     -Features
#        Find symmetry operations
#        Identify space-group type
#        Wyckoff position assignment
#        Refine crystal structure
#        Search irreducible k-points
#        Find a primitive cell
#     -python api for ASE already exists
#       -can probably hijack/rewrite for Project Suite
#




class PointOperators(DevBase):
    
    point_ops_string = '''

nonhex:
E         x, y, z  
C_2x      x,-y,-z  
C_2y     -x, y,-z
C_2z     -x,-y, z
C_2a      y, x,-z
C_2b     -y,-x,-z
C_2c      z,-y, x
C_2d     -x, z, y
C_2e     -z,-y,-x
C_2f     -x,-z,-y
C_31^+    z, x, y
C_32^+   -z, x,-y
C_33^+   -z, x, y
C_34^+    z,-x,-y
C_31^-    y, z, x
C_32^-    y,-z,-x
C_33^-   -y, z,-x
C_34^-   -y,-z, x
C_4x^+    x,-z, y
C_4y^+    z, y,-x
C_4z^+   -y, x, z
C_4x^-    x, z,-y
C_4y^-   -z, y, x
C_4z^-    y,-x, z
I        -x,-y,-z
s_x      -x, y, z
s_y       x,-y, z
s_z       x, y,-z
s_da     -y,-x, z
s_db      y, x, z
s_dc     -z, y,-x
s_dd      x,-z,-y
s_de      z, y, x
s_df      x, z, y
S_61^+   -y,-z,-x
S_62^+   -y, z, x
S_63^+    y,-z, x
S_64^+    y, z,-x
S_61^-   -z,-x,-y
S_62^-    z,-x, y
S_63^-    z, x,-y
S_64^-   -z, x, y
S_4x^+   -x,-z, y
S_4y^+    z,-y,-x
S_4z^+   -y, x,-z
S_4x^-   -x, z,-y
S_4y^-   -z,-y, x
S_4z^-    y,-x,-z                  

hex:
E         x, y, z
C_2      -x,-y, z
C_3^+    -y,x-y,z
C_3^-    y-x,-x,z
C_6^+    x-y,x, z
C_6^-     y,y-x,z
C_21^'   y-x,y,-z
C_22^'    x,x-y,-z
C_23^'   -y,-x,-z
C_21^''  x-y,-y,-z
C_22^''  -x,y-x,-z
C_23^''   y, x,-z
I        -x,-y,-z
s_h       x, y,-z
S_3^+    -y,x-y,-z
S_3^-    y-x,-x,-z
S_6^+    x-y, x,-z
S_6^-     y,y-x,-z
s_d1     x-y,-y,z
s_d2     -x,y-x,z
s_d3      y, x, z
s_v1     y-x,y, z
s_v2      x,x-y,z
s_v3     -y,-x, z

perm:
P_31       x, y, z
P_32       z, x, y
P_33       y, z, x
P_34       x, z, y
P_35       z, y, x
P_36       y, x, z

ref:
R_x       -x, y, z
R_y        x,-y, z
R_z        x, y,-z
R_xy      -x,-y, z
R_xz      -x, y,-z
R_yz       x,-y,-z
R_xyz     -x,-y,-z
'''

#end class PointOperators



class SpaceGroups(DevBase):
    
    space_groups_string = '''

Triclinic:
1         1          E
2        -1          E,I
                     
Monoclinic:          
3-5       2          E,C2z
6-9       m          E,sz
10-15     2/m        E,C2z,I,sz
                     
Orthorhombic:        
16-24     222        E,C2x,C2y,C2z
25-46     mm2        E,C2z,sx,sy
47-74     mmm        E,C2x,C2y,C2z,I,sx,sy,sz
                     
Tetragonal:          
75-80     4          E,C4zp,C4zm,C2z
81-82    -4          E,S4zp,S4zm,C2z
83-88     4/m        E,C4zp,C4zm,C2z,I,S4zp,S4zm,sz
89-98     422        E,C4zp,C4zm,C2z,C2x,C2y,C2a,C2b
99-110    4mm        E,C4zp,C4zm,C2z,sx,sy,sda,sdb
111-122  -42m    (1) E,S4zp,S4zm,C2z,C2x,C2y,sda,sdb
                 (2) E,S4zp,S4zm,C2z,C2a,C2b,sx,sy
123-142   4/mmm      E,C4zp,C4zm,C2z,C2x,C2y,C2a,C2b,I,S4zp,S4zm,sz,sx,sy,sda,sdb

Hexagonal:
143-146   3          E,C3p,C3m
147-148  -3          E,C3p,C3m,I,S6p,S6m
149-155   32     (1) E,C3p,C3m,C21p, C22p, C23p
                 (2) E,C3p,C3m,C21pp,C22pp,C23pp
156-161   3m     (1) E,C3p,C3m,sd1,sd2,sd3
                 (2) E,C3p,C3m,sv1,sv2,sv3
162-176  -3m     (1) E,C3p,C3m,C21p, C22p, C23p, I,S6p,S6m,sd1,sd2,sd3
                 (2) E,C3p,C3m,C21pp,C22pp,C23pp,I,S6p,S6m,sv1,sv2,sv3
168-173   6          E,C6p,C6m,C3p,C3m,C2
174      -6          E,S3p,S3m,C3p,C3m,sh
175-176   6/m        E,C6p,C6m,C3p,C3m,C2,I,S3p,S3m,S6p,S6m,sh
177-182   622        E,C6p,C6m,C3p,C3m,C2,C21p,C22p,C23p,C21pp,C22pp,C23pp
183-186   6mm        E,C6p,C6m,C3p,C3m,C2,sd1,sd2,sd3,sv1,sv2,sv3
187-190  -6m2    (1) E,C3p,C3m,C21p, C22p, C23p, sh,S3p,S3m,sv1,sv2,sv3
                 (2) E,C3p,C3m,C21pp,C22pp,C23pp,sh,S3p,S3m,sd1,sd2,sd3
191-194   6/mmm      E,C6p,C6m,C3p,C3m,C2,C21p,C22p,C23p,C21pp,C22pp,C23pp,
                       I,S3p,S3m,S6p,S6m,sh,sd1,sd2,sd3,sv1,sv2,sv3

Cubic:
195-199   23         E,C31p,C31m,C32p,C32m,C33p,C33m,C34p,C34m,C2x,C2y,C2z
200-206   m-3        E,C31p,C31m,C32p,C32m,C33p,C33m,C34p,C34m,C2x,C2y,C2z,
                       I,S61p,S61m,S62p,S62m,S63p,S63m,S64p,S64m,sx,sy,sz
207-214   432        E,C31p,C31m,C32p,C32m,C33p,C33m,C34p,C34m,C2x,C2y,C2z,
                       C2a,C2b,C2c,C2d,C2e,C2f,C4xp,C4xm,C4yp,C4ym,C4zp,C4zm
215-222  -43m        E,C31p,C31m,C32p,C32m,C33p,C33m,C34p,C34m,C2x,C2y,C2z,
                       sda,sdb,sdc,sdd,sde,sdf,S4xp,S4xm,S4yp,S4ym,S4zp,S4zm
221-230   m-3m       E,C31p,C31m,C32p,C32m,C33p,C33m,C34p,C34m,C2x,C2y,C2z,
                       C2a,C2b,C2c,C2d,C2e,C2f,C4xp,C4xm,C4yp,C4ym,C4zp,C4zm,
                       I,S61p,S61m,S62p,S62m,S63p,S63m,S64p,S64m,sx,sy,sz,
                       sda,sdb,sdc,sdd,sde,sdf,S4xp,S4xm,S4yp,S4ym,S4zp,S4zm
'''

#end class SpaceGroups
