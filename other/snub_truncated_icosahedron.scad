// shape height
sh = 50;
// border width
bw = 2;

include <../polyhedra.scad>;

poly_wire(
  snub_truncated_icosahedron
);

// Elements of the snub_truncated_icosahedron array
//    0 - vertices coordinates
//    1 - face descriptions clockwise
//    2 - vertex kinds
//    3 - face kinds
snub_truncated_icosahedron = [[
  [-0.10336825172505137, -0.08809904942445501, -0.990733850247287], // 0 A vertex
  [-0.15211333253129677, 0.17847739959231437, -0.9721148861636614], // 1 C vertex
  [-0.3921507504891127, 0.2665764490167694, -0.88042875107554], // 2 A vertex
  [-0.5834430876406832, 0.08809904942445503, -0.8073615800710443], // 3 C vertex
  [-0.5346980068344377, -0.17847739959231437, -0.8259805441546699], // 4 A vertex
  [-0.2946605888766219, -0.2665764490167694, -0.9176666792427912], // 5 C vertex
  [0.10336825172505137, 0.08809904942445501, -0.990733850247287], // 6 A vertex
  [0.15211333253129677, -0.17847739959231437, -0.9721148861636614], // 7 C vertex
  [0.3921507504891127, -0.2665764490167694, -0.88042875107554], // 8 A vertex
  [0.5834430876406832, -0.08809904942445503, -0.8073615800710443], // 9 C vertex
  [0.5346980068344377, 0.17847739959231437, -0.8259805441546699], // 10 A vertex
  [0.2946605888766219, 0.2665764490167694, -0.9176666792427912], // 11 C vertex
  [-0.45603584341284925, -0.4338297936655575, -0.7770604993504886], // 12 A vertex
  [-0.6774542972870394, -0.3342235916021083, -0.6552482475397476], // 13 C vertex
  [-0.7770604993504886, -0.4560358434128493, -0.43382979366555746], // 14 A vertex
  [-0.6552482475397476, -0.6774542972870397, -0.3342235916021081], // 15 C vertex
  [-0.43382979366555746, -0.7770604993504888, -0.45603584341284925], // 16 A vertex
  [-0.3342235916021081, -0.6552482475397476, -0.6774542972870394], // 17 C vertex
  [-0.5985830997581744, 0.34573074424110234, -0.7226122924296188], // 18 A vertex
  [-0.3886717985229783, 0.5127009911944226, -0.7655533467114946], // 19 C vertex
  [-0.3457307442411024, 0.7226122924296187, -0.5985830997581744], // 20 A vertex
  [-0.5127009911944226, 0.7655533467114946, -0.3886717985229783], // 21 C vertex
  [-0.7226122924296187, 0.5985830997581743, -0.3457307442411024], // 22 A vertex
  [-0.7655533467114946, 0.38867179852297823, -0.5127009911944226], // 23 C vertex
  [0.45603584341284925, 0.4338297936655575, -0.7770604993504886], // 24 A vertex
  [0.6774542972870394, 0.3342235916021083, -0.6552482475397476], // 25 C vertex
  [0.7770604993504886, 0.4560358434128493, -0.43382979366555746], // 26 A vertex
  [0.6552482475397476, 0.6774542972870397, -0.3342235916021081], // 27 C vertex
  [0.43382979366555746, 0.7770604993504888, -0.45603584341284925], // 28 A vertex
  [0.3342235916021081, 0.6552482475397476, -0.6774542972870394], // 29 C vertex
  [0.5985830997581744, -0.34573074424110234, -0.7226122924296188], // 30 A vertex
  [0.3886717985229783, -0.5127009911944226, -0.7655533467114946], // 31 C vertex
  [0.3457307442411024, -0.7226122924296187, -0.5985830997581744], // 32 A vertex
  [0.5127009911944226, -0.7655533467114946, -0.3886717985229783], // 33 C vertex
  [0.7226122924296187, -0.5985830997581743, -0.3457307442411024], // 34 A vertex
  [0.7655533467114946, -0.38867179852297823, -0.5127009911944226], // 35 C vertex
  [-0.2665764490167694, 0.8804287510755401, -0.3921507504891125], // 36 A vertex
  [-0.08809904942445498, 0.8073615800710445, -0.5834430876406831], // 37 C vertex
  [0.17847739959231437, 0.8259805441546701, -0.5346980068344377], // 38 A vertex
  [0.2665764490167694, 0.9176666792427913, -0.2946605888766218], // 39 C vertex
  [0.08809904942445508, 0.990733850247287, -0.1033682517250513], // 40 A vertex
  [-0.17847739959231435, 0.9721148861636614, -0.1521133325312967], // 41 C vertex
  [0.2665764490167694, -0.8804287510755401, -0.3921507504891125], // 42 A vertex
  [0.08809904942445498, -0.8073615800710445, -0.5834430876406831], // 43 C vertex
  [-0.17847739959231437, -0.8259805441546701, -0.5346980068344377], // 44 A vertex
  [-0.2665764490167694, -0.9176666792427913, -0.2946605888766218], // 45 C vertex
  [-0.08809904942445508, -0.990733850247287, -0.1033682517250513], // 46 A vertex
  [0.17847739959231435, -0.9721148861636614, -0.1521133325312967], // 47 C vertex
  [-0.8259805441546699, -0.5346980068344377, -0.17847739959231446], // 48 A vertex
  [-0.9176666792427912, -0.2946605888766219, -0.26657644901676947], // 49 C vertex
  [-0.990733850247287, -0.10336825172505133, -0.08809904942445514], // 50 A vertex
  [-0.9721148861636614, -0.15211333253129677, 0.17847739959231426], // 51 C vertex
  [-0.88042875107554, -0.3921507504891127, 0.2665764490167693], // 52 A vertex
  [-0.8073615800710443, -0.5834430876406831, 0.08809904942445491], // 53 C vertex
  [-0.990733850247287, 0.10336825172505137, 0.08809904942445501], // 54 A vertex
  [-0.9721148861636614, 0.15211333253129677, -0.17847739959231443], // 55 C vertex
  [-0.88042875107554, 0.3921507504891127, -0.2665764490167694], // 56 A vertex
  [-0.8073615800710443, 0.5834430876406832, -0.08809904942445503], // 57 C vertex
  [-0.8259805441546699, 0.5346980068344377, 0.17847739959231435], // 58 A vertex
  [-0.9176666792427912, 0.2946605888766219, 0.2665764490167694], // 59 C vertex
  [0.8804287510755401, -0.3921507504891125, -0.2665764490167694], // 60 A vertex
  [0.8073615800710445, -0.5834430876406831, -0.08809904942445498], // 61 C vertex
  [0.8259805441546701, -0.5346980068344377, 0.17847739959231437], // 62 A vertex
  [0.9176666792427913, -0.2946605888766218, 0.2665764490167694], // 63 C vertex
  [0.990733850247287, -0.1033682517250513, 0.08809904942445505], // 64 A vertex
  [0.9721148861636614, -0.1521133325312967, -0.17847739959231437], // 65 C vertex
  [0.8259805441546699, 0.5346980068344377, -0.17847739959231446], // 66 A vertex
  [0.9176666792427912, 0.2946605888766219, -0.26657644901676947], // 67 C vertex
  [0.990733850247287, 0.10336825172505133, -0.08809904942445514], // 68 A vertex
  [0.9721148861636614, 0.15211333253129677, 0.17847739959231426], // 69 C vertex
  [0.88042875107554, 0.3921507504891127, 0.2665764490167693], // 70 A vertex
  [0.8073615800710443, 0.5834430876406831, 0.08809904942445491], // 71 C vertex
  [0.08809904942445498, -0.990733850247287, 0.10336825172505137], // 72 A vertex
  [-0.17847739959231443, -0.9721148861636614, 0.15211333253129677], // 73 C vertex
  [-0.2665764490167694, -0.88042875107554, 0.3921507504891127], // 74 A vertex
  [-0.08809904942445503, -0.8073615800710443, 0.5834430876406832], // 75 C vertex
  [0.17847739959231435, -0.8259805441546699, 0.5346980068344377], // 76 A vertex
  [0.2665764490167694, -0.9176666792427912, 0.2946605888766219], // 77 C vertex
  [-0.08809904942445498, 0.990733850247287, 0.10336825172505137], // 78 A vertex
  [0.17847739959231443, 0.9721148861636614, 0.15211333253129677], // 79 C vertex
  [0.2665764490167694, 0.88042875107554, 0.3921507504891127], // 80 A vertex
  [0.08809904942445503, 0.8073615800710443, 0.5834430876406832], // 81 C vertex
  [-0.17847739959231435, 0.8259805441546699, 0.5346980068344377], // 82 A vertex
  [-0.2665764490167694, 0.9176666792427912, 0.2946605888766219], // 83 C vertex
  [-0.34573074424110234, -0.7226122924296188, 0.5985830997581744], // 84 A vertex
  [-0.5127009911944226, -0.7655533467114946, 0.3886717985229783], // 85 C vertex
  [-0.7226122924296187, -0.5985830997581744, 0.3457307442411024], // 86 A vertex
  [-0.7655533467114946, -0.3886717985229783, 0.5127009911944226], // 87 C vertex
  [-0.5985830997581743, -0.34573074424110245, 0.7226122924296188], // 88 A vertex
  [-0.38867179852297823, -0.5127009911944226, 0.7655533467114946], // 89 C vertex
  [-0.7770604993504886, 0.45603584341284925, 0.4338297936655575], // 90 A vertex
  [-0.6552482475397476, 0.6774542972870394, 0.3342235916021083], // 91 C vertex
  [-0.43382979366555746, 0.7770604993504886, 0.4560358434128493], // 92 A vertex
  [-0.3342235916021081, 0.6552482475397476, 0.6774542972870397], // 93 C vertex
  [-0.45603584341284925, 0.43382979366555746, 0.7770604993504886], // 94 A vertex
  [-0.6774542972870394, 0.3342235916021081, 0.6552482475397476], // 95 C vertex
  [0.34573074424110234, 0.7226122924296188, 0.5985830997581744], // 96 A vertex
  [0.5127009911944226, 0.7655533467114946, 0.3886717985229783], // 97 C vertex
  [0.7226122924296187, 0.5985830997581744, 0.3457307442411024], // 98 A vertex
  [0.7655533467114946, 0.3886717985229783, 0.5127009911944226], // 99 C vertex
  [0.5985830997581743, 0.34573074424110245, 0.7226122924296188], // 100 A vertex
  [0.38867179852297823, 0.5127009911944226, 0.7655533467114946], // 101 C vertex
  [0.7770604993504886, -0.45603584341284925, 0.4338297936655575], // 102 A vertex
  [0.6552482475397476, -0.6774542972870394, 0.3342235916021083], // 103 C vertex
  [0.43382979366555746, -0.7770604993504886, 0.4560358434128493], // 104 A vertex
  [0.3342235916021081, -0.6552482475397476, 0.6774542972870397], // 105 C vertex
  [0.45603584341284925, -0.43382979366555746, 0.7770604993504886], // 106 A vertex
  [0.6774542972870394, -0.3342235916021081, 0.6552482475397476], // 107 C vertex
  [-0.3921507504891125, -0.2665764490167694, 0.8804287510755401], // 108 A vertex
  [-0.5834430876406831, -0.08809904942445498, 0.8073615800710445], // 109 C vertex
  [-0.5346980068344377, 0.17847739959231437, 0.8259805441546701], // 110 A vertex
  [-0.2946605888766218, 0.2665764490167694, 0.9176666792427913], // 111 C vertex
  [-0.1033682517250513, 0.08809904942445508, 0.990733850247287], // 112 A vertex
  [-0.1521133325312967, -0.17847739959231435, 0.9721148861636614], // 113 C vertex
  [0.3921507504891125, 0.2665764490167694, 0.8804287510755401], // 114 A vertex
  [0.5834430876406831, 0.08809904942445498, 0.8073615800710445], // 115 C vertex
  [0.5346980068344377, -0.17847739959231437, 0.8259805441546701], // 116 A vertex
  [0.2946605888766218, -0.2665764490167694, 0.9176666792427913], // 117 C vertex
  [0.1033682517250513, -0.08809904942445508, 0.990733850247287], // 118 A vertex
  [0.1521133325312967, 0.17847739959231435, 0.9721148861636614], // 119 C vertex
  [0.2019485620928121, -0.4143275893866615, -0.8839810559279121], // 120 B vertex
  [-0.0443869292354786, -0.3169491449578971, -0.9441642443564814], // 121 B vertex
  [-0.22938119301657425, -0.48614714954590216, -0.8395941266924335], // 122 B vertex
  [-0.09737844442876437, -0.6880957116387143, -0.7147830513399072], // 123 B vertex
  [0.16919800458800502, -0.6437087824032357, -0.7422156822636692], // 124 B vertex
  [-0.2019485620928121, 0.4143275893866615, -0.8839810559279121], // 125 B vertex
  [0.0443869292354786, 0.3169491449578971, -0.9441642443564814], // 126 B vertex
  [0.22938119301657425, 0.48614714954590216, -0.8395941266924335], // 127 B vertex
  [0.09737844442876437, 0.6880957116387143, -0.7147830513399072], // 128 B vertex
  [-0.16919800458800502, 0.6437087824032357, -0.7422156822636692], // 129 B vertex
  [-0.7147830513399072, -0.09737844442876427, -0.6880957116387144], // 130 B vertex
  [-0.7422156822636692, 0.16919800458800507, -0.6437087824032357], // 131 B vertex
  [-0.8839810559279121, 0.2019485620928121, -0.41432758938666153], // 132 B vertex
  [-0.9441642443564814, -0.04438692923547852, -0.31694914495789717], // 133 B vertex
  [-0.8395941266924335, -0.22938119301657425, -0.4861471495459022], // 134 B vertex
  [0.7147830513399072, 0.09737844442876427, -0.6880957116387144], // 135 B vertex
  [0.7422156822636692, -0.16919800458800507, -0.6437087824032357], // 136 B vertex
  [0.8839810559279121, -0.2019485620928121, -0.41432758938666153], // 137 B vertex
  [0.9441642443564814, 0.04438692923547852, -0.31694914495789717], // 138 B vertex
  [0.8395941266924335, 0.22938119301657425, -0.4861471495459022], // 139 B vertex
  [-0.4861471495459022, -0.8395941266924335, -0.22938119301657425], // 140 B vertex
  [-0.6880957116387143, -0.7147830513399072, -0.09737844442876432], // 141 B vertex
  [-0.6437087824032357, -0.7422156822636692, 0.16919800458800505], // 142 B vertex
  [-0.41432758938666153, -0.8839810559279121, 0.20194856209281214], // 143 B vertex
  [-0.31694914495789717, -0.9441642443564814, -0.04438692923547856], // 144 B vertex
  [-0.6437087824032359, 0.7422156822636692, -0.16919800458800502], // 145 B vertex
  [-0.4143275893866616, 0.8839810559279121, -0.20194856209281214], // 146 B vertex
  [-0.3169491449578973, 0.9441642443564814, 0.04438692923547859], // 147 B vertex
  [-0.4861471495459023, 0.8395941266924335, 0.22938119301657428], // 148 B vertex
  [-0.6880957116387144, 0.7147830513399072, 0.09737844442876435], // 149 B vertex
  [0.6437087824032359, -0.7422156822636692, -0.16919800458800502], // 150 B vertex
  [0.4143275893866616, -0.8839810559279121, -0.20194856209281214], // 151 B vertex
  [0.3169491449578973, -0.9441642443564814, 0.04438692923547859], // 152 B vertex
  [0.4861471495459023, -0.8395941266924335, 0.22938119301657428], // 153 B vertex
  [0.6880957116387144, -0.7147830513399072, 0.09737844442876435], // 154 B vertex
  [0.4861471495459022, 0.8395941266924335, -0.22938119301657425], // 155 B vertex
  [0.6880957116387143, 0.7147830513399072, -0.09737844442876432], // 156 B vertex
  [0.6437087824032357, 0.7422156822636692, 0.16919800458800505], // 157 B vertex
  [0.41432758938666153, 0.8839810559279121, 0.20194856209281214], // 158 B vertex
  [0.31694914495789717, 0.9441642443564814, -0.04438692923547856], // 159 B vertex
  [-0.8839810559279121, -0.2019485620928121, 0.4143275893866615], // 160 B vertex
  [-0.9441642443564814, 0.0443869292354786, 0.3169491449578971], // 161 B vertex
  [-0.8395941266924335, 0.22938119301657425, 0.48614714954590216], // 162 B vertex
  [-0.7147830513399072, 0.09737844442876437, 0.6880957116387143], // 163 B vertex
  [-0.7422156822636692, -0.16919800458800502, 0.6437087824032357], // 164 B vertex
  [0.9441642443564814, -0.04438692923547856, 0.31694914495789717], // 165 B vertex
  [0.8395941266924335, -0.22938119301657425, 0.4861471495459022], // 166 B vertex
  [0.7147830513399072, -0.09737844442876432, 0.6880957116387143], // 167 B vertex
  [0.7422156822636692, 0.16919800458800505, 0.6437087824032357], // 168 B vertex
  [0.8839810559279121, 0.20194856209281214, 0.41432758938666153], // 169 B vertex
  [0.09737844442876427, -0.6880957116387144, 0.7147830513399072], // 170 B vertex
  [-0.16919800458800507, -0.6437087824032357, 0.7422156822636692], // 171 B vertex
  [-0.2019485620928121, -0.41432758938666153, 0.8839810559279121], // 172 B vertex
  [0.044386929235478524, -0.31694914495789717, 0.9441642443564814], // 173 B vertex
  [0.22938119301657425, -0.4861471495459022, 0.8395941266924335], // 174 B vertex
  [-0.09737844442876427, 0.6880957116387144, 0.7147830513399072], // 175 B vertex
  [0.16919800458800507, 0.6437087824032357, 0.7422156822636692], // 176 B vertex
  [0.2019485620928121, 0.41432758938666153, 0.8839810559279121], // 177 B vertex
  [-0.044386929235478524, 0.31694914495789717, 0.9441642443564814], // 178 B vertex
  [-0.22938119301657425, 0.4861471495459022, 0.8395941266924335] // 179 B vertex
], [
 [0, 1, 2, 3, 4, 5], // 0 α face
 [6, 7, 8, 9, 10, 11], // 1 α face
 [12, 13, 14, 15, 16, 17], // 2 α face
 [18, 19, 20, 21, 22, 23], // 3 α face
 [24, 25, 26, 27, 28, 29], // 4 α face
 [30, 31, 32, 33, 34, 35], // 5 α face
 [36, 37, 38, 39, 40, 41], // 6 α face
 [42, 43, 44, 45, 46, 47], // 7 α face
 [48, 49, 50, 51, 52, 53], // 8 α face
 [54, 55, 56, 57, 58, 59], // 9 α face
 [60, 61, 62, 63, 64, 65], // 10 α face
 [66, 67, 68, 69, 70, 71], // 11 α face
 [72, 73, 74, 75, 76, 77], // 12 α face
 [78, 79, 80, 81, 82, 83], // 13 α face
 [84, 85, 86, 87, 88, 89], // 14 α face
 [90, 91, 92, 93, 94, 95], // 15 α face
 [96, 97, 98, 99, 100, 101], // 16 α face
 [102, 103, 104, 105, 106, 107], // 17 α face
 [108, 109, 110, 111, 112, 113], // 18 α face
 [114, 115, 116, 117, 118, 119], // 19 α face
 [120, 121, 122, 123, 124], // 20 β face
 [125, 126, 127, 128, 129], // 21 β face
 [130, 131, 132, 133, 134], // 22 β face
 [135, 136, 137, 138, 139], // 23 β face
 [140, 141, 142, 143, 144], // 24 β face
 [145, 146, 147, 148, 149], // 25 β face
 [150, 151, 152, 153, 154], // 26 β face
 [155, 156, 157, 158, 159], // 27 β face
 [160, 161, 162, 163, 164], // 28 β face
 [165, 166, 167, 168, 169], // 29 β face
 [170, 171, 172, 173, 174], // 30 β face
 [175, 176, 177, 178, 179], // 31 β face
 [0, 121, 7], // 32 γ face
 [6, 126, 1], // 33 γ face
 [125, 19, 2], // 34 γ face
 [18, 131, 3], // 35 γ face
 [120, 31, 8], // 36 γ face
 [30, 136, 9], // 37 γ face
 [12, 122, 5], // 38 γ face
 [130, 13, 4], // 39 γ face
 [134, 49, 14], // 40 γ face
 [48, 141, 15], // 41 γ face
 [129, 37, 20], // 42 γ face
 [36, 146, 21], // 43 γ face
 [24, 127, 11], // 44 γ face
 [135, 25, 10], // 45 γ face
 [139, 67, 26], // 46 γ face
 [66, 156, 27], // 47 γ face
 [124, 43, 32], // 48 γ face
 [42, 151, 33], // 49 γ face
 [128, 29, 38], // 50 γ face
 [155, 39, 28], // 51 γ face
 [123, 17, 44], // 52 γ face
 [140, 45, 16], // 53 γ face
 [144, 73, 46], // 54 γ face
 [72, 152, 47], // 55 γ face
 [133, 55, 50], // 56 γ face
 [54, 161, 51], // 57 γ face
 [132, 23, 56], // 58 γ face
 [145, 57, 22], // 59 γ face
 [149, 91, 58], // 60 γ face
 [90, 162, 59], // 61 γ face
 [60, 137, 35], // 62 γ face
 [150, 61, 34], // 63 γ face
 [154, 103, 62], // 64 γ face
 [102, 166, 63], // 65 γ face
 [138, 65, 68], // 66 γ face
 [165, 69, 64], // 67 γ face
 [143, 85, 74], // 68 γ face
 [84, 171, 75], // 69 γ face
 [78, 147, 41], // 70 γ face
 [159, 79, 40], // 71 γ face
 [158, 97, 80], // 72 γ face
 [96, 176, 81], // 73 γ face
 [88, 164, 109], // 74 γ face
 [108, 172, 89], // 75 γ face
 [142, 53, 86], // 76 γ face
 [160, 87, 52], // 77 γ face
 [148, 83, 92], // 78 γ face
 [175, 93, 82], // 79 γ face
 [100, 168, 115], // 80 γ face
 [114, 177, 101], // 81 γ face
 [157, 71, 98], // 82 γ face
 [169, 99, 70], // 83 γ face
 [153, 77, 104], // 84 γ face
 [170, 105, 76], // 85 γ face
 [163, 95, 110], // 86 γ face
 [179, 111, 94], // 87 γ face
 [167, 107, 116], // 88 γ face
 [174, 117, 106], // 89 γ face
 [173, 113, 118], // 90 γ face
 [178, 119, 112], // 91 γ face
 [7, 6, 0], // 92 δ face
 [1, 0, 6], // 93 δ face
 [2, 1, 125], // 94 ε face
 [126, 125, 1], // 95 ζ face
 [19, 18, 2], // 96 δ face
 [3, 2, 18], // 97 δ face
 [4, 3, 130], // 98 ε face
 [131, 130, 3], // 99 ζ face
 [8, 7, 120], // 100 ε face
 [121, 120, 7], // 101 ζ face
 [31, 30, 8], // 102 δ face
 [9, 8, 30], // 103 δ face
 [10, 9, 135], // 104 ε face
 [136, 135, 9], // 105 ζ face
 [5, 4, 12], // 106 δ face
 [13, 12, 4], // 107 δ face
 [14, 13, 134], // 108 ε face
 [130, 134, 13], // 109 ζ face
 [49, 48, 14], // 110 δ face
 [15, 14, 48], // 111 δ face
 [16, 15, 140], // 112 ε face
 [141, 140, 15], // 113 ζ face
 [20, 19, 129], // 114 ε face
 [125, 129, 19], // 115 ζ face
 [37, 36, 20], // 116 δ face
 [21, 20, 36], // 117 δ face
 [22, 21, 145], // 118 ε face
 [146, 145, 21], // 119 ζ face
 [11, 10, 24], // 120 δ face
 [25, 24, 10], // 121 δ face
 [26, 25, 139], // 122 ε face
 [135, 139, 25], // 123 ζ face
 [67, 66, 26], // 124 δ face
 [27, 26, 66], // 125 δ face
 [28, 27, 155], // 126 ε face
 [156, 155, 27], // 127 ζ face
 [32, 31, 124], // 128 ε face
 [120, 124, 31], // 129 ζ face
 [43, 42, 32], // 130 δ face
 [33, 32, 42], // 131 δ face
 [34, 33, 150], // 132 ε face
 [151, 150, 33], // 133 ζ face
 [38, 37, 128], // 134 ε face
 [129, 128, 37], // 135 ζ face
 [29, 28, 38], // 136 δ face
 [39, 38, 28], // 137 δ face
 [40, 39, 159], // 138 ε face
 [155, 159, 39], // 139 ζ face
 [44, 43, 123], // 140 ε face
 [124, 123, 43], // 141 ζ face
 [17, 16, 44], // 142 δ face
 [45, 44, 16], // 143 δ face
 [46, 45, 144], // 144 ε face
 [140, 144, 45], // 145 ζ face
 [73, 72, 46], // 146 δ face
 [47, 46, 72], // 147 δ face
 [50, 49, 133], // 148 ε face
 [134, 133, 49], // 149 ζ face
 [55, 54, 50], // 150 δ face
 [51, 50, 54], // 151 δ face
 [52, 51, 160], // 152 ε face
 [161, 160, 51], // 153 ζ face
 [56, 55, 132], // 154 ε face
 [133, 132, 55], // 155 ζ face
 [23, 22, 56], // 156 δ face
 [57, 56, 22], // 157 δ face
 [58, 57, 149], // 158 ε face
 [145, 149, 57], // 159 ζ face
 [91, 90, 58], // 160 δ face
 [59, 58, 90], // 161 δ face
 [35, 34, 60], // 162 δ face
 [61, 60, 34], // 163 δ face
 [62, 61, 154], // 164 ε face
 [150, 154, 61], // 165 ζ face
 [103, 102, 62], // 166 δ face
 [63, 62, 102], // 167 δ face
 [64, 63, 165], // 168 ε face
 [166, 165, 63], // 169 ζ face
 [68, 67, 138], // 170 ε face
 [139, 138, 67], // 171 ζ face
 [65, 64, 68], // 172 δ face
 [69, 68, 64], // 173 δ face
 [70, 69, 169], // 174 ε face
 [165, 169, 69], // 175 ζ face
 [74, 73, 143], // 176 ε face
 [144, 143, 73], // 177 ζ face
 [85, 84, 74], // 178 δ face
 [75, 74, 84], // 179 δ face
 [76, 75, 170], // 180 ε face
 [171, 170, 75], // 181 ζ face
 [41, 40, 78], // 182 δ face
 [79, 78, 40], // 183 δ face
 [80, 79, 158], // 184 ε face
 [159, 158, 79], // 185 ζ face
 [97, 96, 80], // 186 δ face
 [81, 80, 96], // 187 δ face
 [82, 81, 175], // 188 ε face
 [176, 175, 81], // 189 ζ face
 [86, 85, 142], // 190 ε face
 [143, 142, 85], // 191 ζ face
 [53, 52, 86], // 192 δ face
 [87, 86, 52], // 193 δ face
 [109, 108, 88], // 194 δ face
 [89, 88, 108], // 195 δ face
 [92, 91, 148], // 196 ε face
 [149, 148, 91], // 197 ζ face
 [83, 82, 92], // 198 δ face
 [93, 92, 82], // 199 δ face
 [94, 93, 179], // 200 ε face
 [175, 179, 93], // 201 ζ face
 [98, 97, 157], // 202 ε face
 [158, 157, 97], // 203 ζ face
 [71, 70, 98], // 204 δ face
 [99, 98, 70], // 205 δ face
 [115, 114, 100], // 206 δ face
 [101, 100, 114], // 207 δ face
 [104, 103, 153], // 208 ε face
 [154, 153, 103], // 209 ζ face
 [77, 76, 104], // 210 δ face
 [105, 104, 76], // 211 δ face
 [106, 105, 174], // 212 ε face
 [170, 174, 105], // 213 ζ face
 [110, 109, 163], // 214 ε face
 [164, 163, 109], // 215 ζ face
 [95, 94, 110], // 216 δ face
 [111, 110, 94], // 217 δ face
 [112, 111, 178], // 218 ε face
 [179, 178, 111], // 219 ζ face
 [116, 115, 167], // 220 ε face
 [168, 167, 115], // 221 ζ face
 [107, 106, 116], // 222 δ face
 [117, 116, 106], // 223 δ face
 [118, 117, 173], // 224 ε face
 [174, 173, 117], // 225 ζ face
 [113, 112, 118], // 226 δ face
 [119, 118, 112], // 227 δ face
 [0, 5, 121], // 228 ε face
 [122, 121, 5], // 229 ζ face
 [12, 17, 122], // 230 ε face
 [123, 122, 17], // 231 ζ face
 [6, 11, 126], // 232 ε face
 [127, 126, 11], // 233 ζ face
 [24, 29, 127], // 234 ε face
 [128, 127, 29], // 235 ζ face
 [18, 23, 131], // 236 ε face
 [132, 131, 23], // 237 ζ face
 [30, 35, 136], // 238 ε face
 [137, 136, 35], // 239 ζ face
 [60, 65, 137], // 240 ε face
 [138, 137, 65], // 241 ζ face
 [48, 53, 141], // 242 ε face
 [142, 141, 53], // 243 ζ face
 [36, 41, 146], // 244 ε face
 [147, 146, 41], // 245 ζ face
 [78, 83, 147], // 246 ε face
 [148, 147, 83], // 247 ζ face
 [42, 47, 151], // 248 ε face
 [152, 151, 47], // 249 ζ face
 [72, 77, 152], // 250 ε face
 [153, 152, 77], // 251 ζ face
 [66, 71, 156], // 252 ε face
 [157, 156, 71], // 253 ζ face
 [54, 59, 161], // 254 ε face
 [162, 161, 59], // 255 ζ face
 [90, 95, 162], // 256 ε face
 [163, 162, 95], // 257 ζ face
 [88, 87, 164], // 258 ε face
 [160, 164, 87], // 259 ζ face
 [102, 107, 166], // 260 ε face
 [167, 166, 107], // 261 ζ face
 [100, 99, 168], // 262 ε face
 [169, 168, 99], // 263 ζ face
 [84, 89, 171], // 264 ε face
 [172, 171, 89], // 265 ζ face
 [108, 113, 172], // 266 ε face
 [173, 172, 113], // 267 ζ face
 [96, 101, 176], // 268 ε face
 [177, 176, 101], // 269 ζ face
 [114, 119, 177], // 270 ε face
 [178, 177, 119] // 271 ζ face
], [
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
], [
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 5, 3, 3, 4, 5,
  4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5,
  3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5,
  4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5,
  3, 3, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3,
  4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 3, 3, 4, 5, 3, 3,
  4, 5, 4, 5, 3, 3, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 3, 3, 4, 5,
  4, 5, 3, 3, 4, 5, 3, 3, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5,
  4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5,
  4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5
]];
