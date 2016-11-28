package by.vma.lab5;

public class Main {
    private static class Matrix {
        public double[][] matrix;
        private int lines;
        private int columns;

        public Matrix(int lines, int columns) throws Exception {
            if (lines < 1 || columns < 1) {
                throw new Exception("Неверный размер.");
            }
            this.lines = lines;
            this.columns = columns;
            this.matrix = new double[lines][columns];
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

        public void fillE() {
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
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

        public Matrix mul(double num) throws Exception {
            Matrix result = new Matrix(lines, columns);
            for(int i = 0; i < lines; i++) {
                for(int j = 0; j < columns; j++) {
                    result.matrix[i][j] = this.matrix[i][j] * num;
                }
            }
            return result;
        }

        public Matrix subtract(Matrix mtr) throws Exception {
            Matrix result = new Matrix(lines, columns);
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    result.matrix[i][j] = this.matrix[i][j] - mtr.matrix[i][j];
                }
            }
            return result;
        }
    }

    private static class Vector {
        public double[] vector;
        private int length;

        public Vector(int length) throws Exception {
            if (length < 1) {
                throw new Exception("Неверный размер.");
            }
            this.length = length;
            vector = new double[length];
        }

        public Vector(Vector init) throws Exception {
            this(init.getLength());
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

        public Vector mul(double num) throws Exception {
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] * num;
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
    private static Vector b;
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
            x = methodJacobi();
            if (x != null) {
                System.out.println("Вектор X(Метод Якоби): ");
                x.print(false);
                System.out.println();
                r = A.mul(x).subtract(b);
                System.out.println("Вектор невязки R(Метод Якоби): ");
                r.print(true);
                System.out.println();
                System.out.println("Норма вектора невязки ||R||(Метод Якоби) = " + r.normI());
            } else {
                System.out.println("Метод Якоби непременим к данной матрице.");
            }
            System.out.println();
            x = methodSeidel();
            if (x != null) {
                System.out.println("Вектор X(Метод Зейделя): ");
                x.print(false);
                System.out.println();
                r = A.mul(x).subtract(b);
                System.out.println("Вектор невязки R(Метод Зейделя): ");
                r.print(true);
                System.out.println();
                System.out.println("Норма вектора невязки ||R||(Метод Зейделя) = " + r.normI());
            } else {
                System.out.println("Метод Зейделя непременим к данной матрице.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static Vector methodJacobi() throws Exception {
        Vector prevX;
        Vector nextX = new Vector(n);
        Matrix B = new Matrix(n, n);
        int counter = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    B.matrix[i][j] = - A.matrix[i][j] / A.matrix[i][i];
                } else {
                    B.matrix[i][j] = 0;
                }
            }
        }
        if (isMethodJacobiConverges()) {
            for (int i = 0; i < n; i++) {
                nextX.vector[i] = b.vector[i] / A.matrix[i][i];
            }
            do {
                prevX = new Vector(nextX);
                for (int i = 0; i < n; i++) {
                    nextX.vector[i] = b.vector[i] / A.matrix[i][i];
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            nextX.vector[i] -= (A.matrix[i][j] / A.matrix[i][i]) * prevX.vector[j];
                        }
                    }
                }
                counter++;
            } while (nextX.subtract(prevX).normI() > epsilon);
            System.out.println("Априорная оценка метода Якоби: " +
                    ((int)((Math.log(epsilon * (1 - B.normI()) / b.normI())) / Math.log(B.normI())) + 1));
            System.out.println("Количество итераций в методе Якоби: " + counter);
            System.out.println();
            return nextX;
        } else {
            return null;
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
                if (i != j) {
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

    private static Vector methodSeidel() throws Exception {
        double tetta = 1 / A.normI();
        int counter = 0;
        Vector g = b.mul(tetta);
        Vector prevX;
        Vector nextX = new Vector(g);
        Matrix B = new Matrix(n, n);
        B.fillE();
        B = B.subtract(A.mul(tetta));
        if (isMethodSeidelConverges(B)) {
            do {
                prevX = new Vector(nextX);
                for (int i = 0; i < n; i++) {
                    nextX.vector[i] = g.vector[i];
                    for (int j = 0; j < i; j++) {
                        nextX.vector[i] += B.matrix[i][j] * nextX.vector[j];
                    }
                    for(int j = i; j < n; j++) {
                        nextX.vector[i] += B.matrix[i][j] * prevX.vector[j];
                    }
                }
                counter++;
            } while (nextX.subtract(prevX).normI() > epsilon);
            System.out.println("Априорная оценка метода Зейделя: " +
                    ((int)((Math.log(epsilon * (1 - B.normI()) / g.normI())) / Math.log(B.normI())) + 1));
            System.out.println("Количество итераций в методе Зейделя: " + counter);
            System.out.println();
            return nextX;
        } else {
            return null;
        }
    }

    private static boolean isMethodSeidelConverges(Matrix B) {
        double sum;
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                sum += Math.abs(B.matrix[i][j]);
            }
            if (sum < 1) {
                return true;
            }
        }
        for (int j = 0; j < n; j++) {
            sum = 0;
            for (int i = 0; i < n; i++) {
                sum += Math.abs(B.matrix[i][j]);
            }
            if (sum < 1) {
                return true;
            }
        }
        return false;
    }
}
