//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define constants required while performing numerical integration
// Values obtained from the listed references
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_CONSTANTS
#define HEADERFILE_CONSTANTS

// Define pi
const double pi = 3.14159265358979323846;

// Define constants used in the symplectic RK4 method
const double g1 = 0.0, g2 = 0.205177662, g3 = 0.608198943, g4 = 0.487278067, g5 = 1.0;
const double b1 = 0.0617588581, b2 = 0.338978026, b3 = 0.614791307, b4 = -0.140548014, b5 = 0.125019823;
const double a21 = b1 * (g2 - g1);
const double a31 = b1 * (g3 - g1), a32 = b2 * (g3 - g2);
const double a41 = b1 * (g4 - g1), a42 = b2 * (g4 - g2), a43 = b3 * (g4 - g3);
const double a51 = b1 * (g5 - g1), a52 = b2 * (g5 - g2), a53 = b3 * (g5 - g3), a54 = b4 * (g5 - g4);
const double B1 = b1 * (1.0 - g1), B2 = b2 * (1.0 - g2), B3 = b3 * (1.0 - g3), B4 = b4 * (1.0 - g4), B5 = b5 * (1.0 - g5);

// Define constants used in the symplectic RK6 method
const double g64 = 0.5, g65 = 0.06520862987680341024, g66 = 0.65373769483744778901, g67 = 0.05586607811787376572;
const double g61 = 1.0 - g67, g62 = 1.0 - g66, g63 = 1.0 - g65;
const double b61 = -0.68774007118557290171, b62 = 0.13118241020105280626, b63 = 0.92161977504885189358, b64 = 0.26987577187133640373;
const double b65 = b63, b66 = b62, b67 = b61;
const double a621 = b61 * (g62 - g61);
const double a631 = b61 * (g63 - g61), a632 = b62 * (g63 - g62);
const double a641 = b61 * (g64 - g61), a642 = b62 * (g64 - g62), a643 = b63 * (g64 - g63);
const double a651 = b61 * (g65 - g61), a652 = b62 * (g65 - g62), a653 = b63 * (g65 - g63), a654 = b64 * (g65 - g64);
const double a661 = b61 * (g66 - g61), a662 = b62 * (g66 - g62), a663 = b63 * (g66 - g63), a664 = b64 * (g66 - g64), a665 = b65 * (g66 - g65);
const double a671 = b61 * (g67 - g61), a672 = b62 * (g67 - g62), a673 = b63 * (g67 - g63), a674 = b64 * (g67 - g64), a675 = b65 * (g67 - g65), a676 = b66 * (g67 - g66);
const double B61 = b61 * (1.0 - g61), B62 = b62 * (1.0 - g62), B63 = b63 * (1.0 - g63), B64 = b64 * (1.0 - g64), B65 = b65 * (1.0 - g65), B66 = b66 * (1.0 - g66), B67 = b67 * (1.0 - g67);

#endif