package by.vma.lab5;

public class Main {
    private static class Matrix {
        public double[][] matrix;
        private int lines;
        private int columns;

        public Matrix() {
            this.lines = 0;
            this.columns = 0;
            this.matrix = null;
        }

        public Matrix(int lines, int columns) throws Exception {
            if (lines < 1 || columns < 1) {
                throw new Exception("Неверный размер.");
            }
            this.lines = lines;
            this.columns = columns;
            this.matrix = new double[lines][columns];
        }

        public Matrix(Matrix init) throws Exception {
            this(init.getLines(), init.getColumns());
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    this.matrix[i][j] = init.matrix[i][j];
                }
            }
        }

        public int getLines() {
            return lines;
        }

        public int getColumns() {
            return columns;
        }

        public void setLine(int index, Vector line) throws Exception {
            if (line.getLength() != columns) {
                throw new Exception("Неверный вектор.");
            }
            for (int i = 0; i < columns; i++) {
                matrix[index][i] = line.vector[i];
            }
        }

        public void setColumn(int index, Vector column) throws Exception {
            if (column.getLength() != lines) {
                throw new Exception("Неверный вектор.");
            }
            for (int i = 0; i < lines; i++) {
                matrix[i][index] = column.vector[i];
            }
        }

        public void print() {
            for (double[] i : matrix) {
                for (double j : i) {
                    System.out.printf("%.5f", j);
                    System.out.print("  ");
                }
                System.out.println();
            }
        }

        public void swap(int fi, int fj, int si, int sj) {
            double tmp = matrix[fi][fj];
            matrix[fi][fj] = matrix[si][sj];
            matrix[si][sj] = tmp;
        }

        public void swapLines(int fline, int sline) {
            for (int i = 0; i < columns; i++) {
                swap(fline, i, sline, i);
            }
        }

        public void swapColumns(int fcolumn, int scolumn) {
            for (int i = 0; i < lines; i++) {
                swap(i, fcolumn, i, scolumn);
            }
        }

        public void fillDefault() {
            double[][] a = {{0.6444, 0.0000, -0.1683, 0.1184, 0.1973},
                    {-0.0395, 0.4208, 0.0000, -0.0802, 0.0263},
                    {0.0132, -0.1184, 0.7627, 0.0145, 0.0460},
                    {0.0395, 0.0000, -0.0960, 0.7627, 0.0000},
                    {0.0263, -0.0395, 0.1907, -0.0158, 0.5523}};
            this.lines = 5;
            this.columns = 5;
            this.matrix = a;
        }

        public void fillE(int n) {
            lines = n;
            columns = n;
            matrix = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        matrix[i][j] = 1;
                    } else {
                        matrix[i][j] = 0;
                    }
                }
            }
        }

        public double normI() {
            double maxSum = -1;
            double sum;
            for (int i = 0; i < lines; i++) {
                sum = 0;
                for (int j = 0; j < columns; j++) {
                    sum += Math.abs(matrix[i][j]);
                }
                if (sum > maxSum) {
                    maxSum = sum;
                }
            }
            return maxSum;
        }

        public Vector mul(Vector vector) throws Exception {
            if (columns != vector.getLength()) {
                throw new Exception("Неверная матрица или вектор.");
            }
            Vector result = new Vector(vector.getLength());
            for (int i = 0; i < lines; i++) {
                result.vector[i] = 0;
                for (int j = 0; j < columns; j++) {
                    result.vector[i] += matrix[i][j] * vector.vector[j];
                }
            }
            return result;
        }

        public Matrix mul(Matrix mtr) throws Exception {
            if (columns != mtr.getLines()) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(lines, mtr.getColumns());
            for (int i = 0; i < result.getLines(); i++) {
                for (int j = 0; j < result.getColumns(); j++) {
                    result.matrix[i][j] = 0;
                    for (int k = 0; k < columns; k++) {
                        result.matrix[i][j] += this.matrix[i][k] * mtr.matrix[k][j];
                    }
                }
            }
            return result;
        }

        public Matrix subtract(Matrix mtr) throws Exception {
            if (lines != mtr.getLines() || columns != mtr.getColumns()) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(lines, columns);
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    result.matrix[i][j] = this.matrix[i][j] - mtr.matrix[i][j];
                }
            }
            return result;
        }

        public Matrix transpose() throws Exception {
            if (lines != columns) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(this);
            for (int i = 0; i < lines; i++) {
                for (int j = i + 1; j < columns; j++) {
                    result.swap(i, j, j, i);
                }
            }
            return result;
        }
    }

    private static class Vector {
        public double[] vector;
        private int length;

        public Vector() {
            this.length = 0;
            this.vector = null;
        }

        public Vector(int length) throws Exception {
            if (length < 1) {
                throw new Exception("Неверный размер.");
            }
            this.length = length;
            vector = new double[length];
        }

        public Vector(Vector init) {
            this.length = init.getLength();
            for (int i = 0; i < length; i++) {
                this.vector[i] = init.vector[i];
            }
        }

        public int getLength() {
            return length;
        }

        public void print(boolean exponent) {
            for (double item : vector) {
                if (exponent) {
                    System.out.printf("%e\n", item);
                } else {
                    System.out.printf("%.5f\n", item);
                }
            }
        }

        public void swap(int i, int j) {
            double tmp = vector[i];
            vector[i] = vector[j];
            vector[j] = tmp;
        }

        public void fillDefault() {
            double[] b = {1.2677, 1.6819, -2.3657, -6.5369, 2.8351};
            this.length = 5;
            this.vector = b;
        }

        public Vector subtract(Vector sub) throws Exception {
            if (length != sub.getLength()) {
                throw new Exception("Неверный вектор.");
            }
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] - sub.vector[i];
            }
            return result;
        }

        public Vector add(Vector second) throws Exception {
            if (length != second.getLength()) {
                throw new Exception("Неверный вектор.");
            }
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] + second.vector[i];
            }
            return result;
        }

        public double normI() {
            double max = Math.abs(vector[0]);
            for (int i = 1; i < length; i++) {
                if (Math.abs(vector[i]) > max) {
                    max = Math.abs(vector[i]);
                }
            }
            return max;
        }
    }

    private static Matrix A;
    private static Matrix B;
    private static Vector b;
    private static Vector g;
    private static int n = 5;
    private static double epsilon = 0.00001;

    public static void main(String[] args) {
        Vector x;
        Vector r;
        try {
            A = new Matrix(n, n);
            b = new Vector(n);
            A.fillDefault();
            b.fillDefault();
            if (isMethodJacobiConverges()) {
                x = methodJacobi();
                System.out.println("Вектор X: ");
                x.print(false);
                System.out.println();
                r = A.mul(x).subtract(b);
                System.out.println("Вектор невязки R: ");
                r.print(true);
                System.out.println();
            } else {
                System.out.println("Метод Якоби непременим к данной матрице.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static Vector methodJacobi() throws Exception {
        Vector prevX;
        Vector nextX;
        fillBandG();
        nextX = g;
        do {
            prevX = nextX;
            nextX = B.mul(prevX).add(g);
        } while (nextX.subtract(prevX).normI() > epsilon);
        return nextX;
    }

    private static void fillBandG() throws Exception {
        B = new Matrix(n, n);
        g = new Vector(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    B.matrix[i][j] = -A.matrix[i][j] / A.matrix[i][i];
                }
            }
            g.vector[i] = b.vector[i] / A.matrix[i][i];
        }
    }

    private static boolean isMethodJacobiConverges() {
        double sum;
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                sum += Math.abs(A.matrix[i][j] / A.matrix[i][i]);
            }
            if (sum < 1) {
                return true;
            }
        }
        for (int j = 0; j < n; j++) {
            sum = 0;
            for (int i = 0; i < n; i++) {
                sum += Math.abs(A.matrix[i][j] / A.matrix[i][i]);
            }
            if (sum < 1) {
                return true;
            }
        }
        sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if(i != j) {
                    sum += Math.pow(A.matrix[i][j] / A.matrix[i][i], 2);
                }
            }
        }
        if (sum < 1) {
            return true;
        } else {
            return false;
        }
    }
}
