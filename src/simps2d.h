#ifndef SIMPS2D_SEEN
#define SIMPS2D_SEEN

// Integration 2D domaine carre
// Selon Abramowitz p892-893 formule 25.462, residu ordre 0.5**6
// 0.387298=0.5*sqrt(3/5)  0.197531=16/81  0.123457=10/81  0.077160=25/324
#if defined(SIMPSON9)
  static int   nd2d = 9;
  static float dx2d[9] = { 0.000000 ,  0.000000 ,  0.000000
                         ,-0.387298 ,  0.387298 , -0.387298
                         , 0.387298 , -0.387298 ,  0.387298 };
  static float dy2d[9] = { 0.000000 , -0.387298 ,  0.387298
                         , 0.000000 ,  0.000000 , -0.387298
                         , 0.387298 ,  0.387298 , -0.387298 };
  static float  w2d[9] = { 0.197531 ,  0.123457 ,  0.123457
                         , 0.123457 ,  0.123457 ,  0.077157
                         , 0.0771605,  0.0771605,  0.0771605};
#elif defined(SIMPSON4)
  static int   nd2d = 4;
  static float dx2d[4] = { 0.288675, 0.288675,-0.288675,-0.288675 };
  static float dy2d[4] = { 0.288675,-0.288675, 0.288675,-0.288675 };
  static float  w2d[4] = { 0.250000, 0.250000, 0.250000, 0.250000 };
#elif defined(INTEG5)
  static int   nd2d = 5;
  static float dx2d[5] = { 0. , -0.3, -0.3,  0.3, 0.3 };
  static float dy2d[5] = { 0. , -0.3,  0.3, -0.3, 0.3 };
  static float  w2d[5] = { 0.2,  0.2,  0.2,  0.2, 0.2 };
#elif  defined(NOINTEG)
  static int   nd2d = 1;
  static float dx2d[1] = { 0. };
  static float dy2d[1] = { 0. };
  static float  w2d[1] = { 1. };
#else
  static int   nd2d = 0;
  static float dx2d[1] = { 999999. };
  static float dy2d[1] = { 999999. };
  static float  w2d[1] = { 0. };
#endif

#endif
