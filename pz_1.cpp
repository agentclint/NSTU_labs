#include <math.h> 
#include <stdio.h> 

int main() 
{ 
double h, y, t = 0; 
h = 0.1; 
y = 1.0; // начальное значение y_n 
for (t = 0.0; t < 1; t += h) { 
printf("%lf\n", y); 
y += h * 2 * t * y; 
} 
return 0; 
}