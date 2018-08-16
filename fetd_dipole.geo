gap_leng = 0.2;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {4, 0, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};
Point(5) = {0, 4, 0, 1.0};
Point(6) = {2, 2, 0, 1.0};
Point(7) = {4, 4, 0, 1.0};
Point(8) = {4, 2, 0, 1.0};
Point(9) = {2, 4, 0, 1.0};
Line(1) = {5, 9};
Line(2) = {9, 7};
Line(3) = {7, 8};
Line(4) = {8, 6};
Line(5) = {6, 9};
Line(6) = {5, 4};
Line(7) = {4, 6};
Line(8) = {6, 2};
Line(9) = {2, 1};
Line(10) = {1, 4};
Line(11) = {2, 2};
Line(12) = {3, 2};
Line(13) = {3, 8};
Line Loop(14) = {6, 7, 5, -1};
Plane Surface(15) = {14};
Line Loop(16) = {5, 2, 3, 4};
Plane Surface(17) = {16};
Line Loop(18) = {7, 8, 9, 10};
Plane Surface(19) = {18};
Line Loop(20) = {4, 8, -12, 13};
Plane Surface(21) = {20};
Extrude {0, 0, 2 - gap_leng / 2} {
  Surface{15, 17, 19, 21};
}
Extrude {0, 0, gap_leng} {
  Surface{43, 65, 87, 109};
}
Extrude {0, 0, 2 - gap_leng / 2} {
  Surface{175, 131, 153, 197};
}
Physical Line(286) = {121};
Physical Volume(287) = {7, 5, 9, 10, 3, 1, 4, 2, 8, 6, 12, 11};

Physical Surface(290) = {60, 258, 284, 108, 196, 148, 104, 82, 15, 214, 280, 192, 170, 218, 228, 86, 30, 118, 174, 240, 254, 56, 42, 130, 144, 241, 263, 285, 219, 21, 19, 17};
