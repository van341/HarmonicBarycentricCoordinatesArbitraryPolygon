
using barycentric_calc;
using System.Drawing;
using System.Numerics;
using System.Drawing.Imaging;
/*
 This is a test example of calculating the external and internal barycentric coordinates for an 
arbitrary polygon using an approximate analytical method developed by Ivan Sergeevich Polyanskii.
https://www.researchgate.net/profile/Ivan-Polyanskii
https://www.mathnet.ru/php/person.phtml?option_lang=rus&personid=117188&ysclid=lz7a0p40l5655517880
*/


//=========configurable parameters=======
// the order of the polynomial of approximate analytical calculation of barycentric coordinates
int NK = 8;
// the number of integration segments on the edge
int MN = 4;
// the number of integration nodes on the segment
int MM = 24;
//the linear size of the rectangular area in which the calculation will be performed
double a = 10;
//the number of points in the calculated grid along one axis
int M = 500;
//the number of the barycentric coordinate whose values will be saved to a file
int k = 6;
//the path to the file containing the coordinates of the vertices of the polygon under study (specify the path to the file!!!)
string FileNameCoordinatesVerticesPolygon = "...\\HarmonicBarycentricCoordinatesArbitraryPolygon\\CoordinatesVerticesPolygon.txt";
//the path and name of the file in which the calculation results will be saved - internal barycentric coordinates (specify the path to the file!!!)
string FileNameResultImageIBC = "...\\HarmonicBarycentricCoordinatesArbitraryPolygon\\ResultImageEBC.png";
//the path and name of the file in which the calculation results will be saved  - external barycentric coordinates (specify the path to the file!!!)
string FileNameResultImageEBC = "...\\HarmonicBarycentricCoordinatesArbitraryPolygon\\ResultImageIBC.png";

//==============calculations============
//reading the coordinates of the vertices of an arbitrary polygon from a file
static List<Complex> ReadDoubleArrayFromFile(string filePath)
{
    var Point = new List<Complex>();
    foreach (var line in File.ReadLines(filePath))
    {
        var elements = line.Split('\t');
        Point.Add(new Complex(double.Parse(elements[0]), double.Parse(elements[1])));
    }
    return Point;
}

//the procedure for saving the calculation results to a file
static void SaveResultToFile(double[][] result, string filePath)
{
    using (StreamWriter writer = new StreamWriter(filePath))
    {
        foreach (var array in result)
        {
            writer.WriteLine(string.Join("\t", array));
        }
    }
}

//the procedure for saving the calculation results to an image
static void SaveResultToImage(double[][] result, int M, string filePath)
{
    using (Bitmap bitmap = new Bitmap(M, M))
    {
        Parallel.For(0, M, i =>
        {
            for (int j = 0; j < M; j++)
            {
                double zetaValue = result[i * M + j][2];
                double normalizedValue = (zetaValue + 0.1) / 1.2;
                Color color = ColorFromHSV(normalizedValue * 360, 1.0, 1.0);
                lock (bitmap)
                {
                    bitmap.SetPixel(i, j, color);
                }
            }
        });
        bitmap.Save(filePath, ImageFormat.Png);
    }
}

//rainbow graph color maps
static Color ColorFromHSV(double hue, double saturation, double value)
{
    hue = (hue + 180) % 360;
    int hi = Convert.ToInt32(Math.Floor(hue / 60)) % 6;
    double f = hue / 60 - Math.Floor(hue / 60);

    value = value * 255;
    int v = Convert.ToInt32(value);
    int p = Convert.ToInt32(value * (1 - saturation));
    int q = Convert.ToInt32(value * (1 - f * saturation));
    int t = Convert.ToInt32(value * (1 - (1 - f) * saturation));

    return hi switch
    {
        0 => Color.FromArgb(255, v, t, p),
        1 => Color.FromArgb(255, q, v, p),
        2 => Color.FromArgb(255, p, v, t),
        3 => Color.FromArgb(255, p, q, v),
        4 => Color.FromArgb(255, t, p, v),
        _ => Color.FromArgb(255, v, p, q),
    };
}

//the vertices of the polygon, set in the order of positive traversal (counterclockwise)
var Point = ReadDoubleArrayFromFile(FileNameCoordinatesVerticesPolygon);
//creating an object for calculating external and internal barycentric coordinates for an arbitrary polygon
var BC_int_ext = new barycentric_generalized_calculations(Point, NK, MN, MM);
//calculation of barycentric coordinates in a given area
double[][] resultEBC = new double[M * M][];
double[][] resultIBC = new double[M * M][];
Parallel.For(0, M, i =>
{
    double x = -0.5 * a + a * i / (M - 1);
    for (int j = 0; j < M; j++)
    {
        double y = -0.5 * a + a * j / (M - 1);
        var P = new Complex(x, y);
        bool checkP = BC_int_ext.check_point_in_polygon(P);
        List<double> zeta = checkP ? BC_int_ext.barycentric_interier(P) : BC_int_ext.barycentric_exterier(P);
        double EBC = checkP ? zeta[k] : -0.1;
        double IBC = checkP ? -0.1 : zeta[k];
        resultEBC[i * M + j] = new double[] { x, y, EBC };
        resultIBC[i * M + j] = new double[] { x, y, IBC };
    }
});

//saving calculation results to an image
SaveResultToImage(resultEBC, M, FileNameResultImageEBC);
SaveResultToImage(resultIBC, M, FileNameResultImageIBC);