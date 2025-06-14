package util;

import util.AClass;

import java.util.List;

public class CombinationChecker {
    private List<AClass> classes;
    private List<String> courses;
    int nbCourses;
    int nbClasses;
    private SolutionClass[] x; //x[i] is the solution for class i
    private boolean[][] overlap; // overlap[i][j] = true means that class i and j overlap
    List<Integer>[] classIndicesOfCourse; // classIndicesOfCourse[i] - the list of classes in course i
    int[] courseOfClass;
    private int idxClass; // the current class being checked
    private int idxCourse; // the course of the class being checked
    private int[] y; // the current combination
    private boolean ans; // whether or not the class being checked belongs to a combination

    private boolean check(int v, int i) { // objective Y[0], Y[1], .. pair-wise not-overlap
        //check if class v overlaps with any of the classes put in the current combination
        for (int j = 0; j <= i - 1; j++) {
            if (overlap[y[j]][v]) return false;
        }
        return true;
    }

    private void solution() {
        ans = true;
    }

    private void tryY(int i) {// try all values for Y[i]
        if (ans) return;
        if (i == idxCourse) {
            if (check(idxClass, i)) {
                y[i] = idxClass;
                if (i == nbCourses - 1) {
                    solution();
                } else {
                    tryY(i + 1);
                }
            }
            return;
        }
        for (int v : classIndicesOfCourse[i]) {
            if (check(v, i)) {
                y[i] = v;
                if (i == nbCourses - 1) {
                    solution();
                } else {
                    tryY(i + 1);
                }
            }
        }
    }

    public CombinationChecker(List<Integer>[] classIndicesOfCourse, List<AClass> classes, SolutionClass[] x) {
        this.classIndicesOfCourse = classIndicesOfCourse;
        this.classes = classes;
        this.x = x;
        nbCourses = classIndicesOfCourse.length;
        nbClasses = classes.size();
        courseOfClass = new int[nbClasses];
        for (int i = 0; i < nbCourses; i++) {
            for (int j : classIndicesOfCourse[i]) {
                courseOfClass[j] = i;
            }
        }
        overlap = new boolean[classes.size()][classes.size()];
        for (int i = 0; i < classes.size(); i++) {
            for (int j = 0; j < classes.size(); j++) {
                overlap[i][j] = x[i].overlap(x[j]);
            }
        }
    }

    public boolean checkInCombination(int idxClass) {
        this.idxClass = idxClass;
        idxCourse = courseOfClass[idxClass];
        ans = false;
        y = new int[nbCourses];// Y[i] index class of course i selected in the combination
        tryY(0);
        return ans;
    }
}
