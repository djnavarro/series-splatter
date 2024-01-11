// variation on the "skipping stone" cellular automaton

#include <Rcpp.h>
using namespace Rcpp;

// wrap position to grid
int wrap(int pos, int size) {
  if(pos < 0) pos = pos + size;
  if(pos >= size) pos = pos - size;
  return pos;
}

// automaton run function
// [[Rcpp::export]]
NumericMatrix automaton(int n_rows, int n_cols, int iterations, int max_span) {

  int source_row = 0;
  int source_col = 0;
  int span_row = 0;
  int span_col = 0;
  int row = 0;
  int col = 0;
  int r = 0;
  int c = 0;
  double source_val = 0;

  NumericMatrix grid(n_rows, n_cols);
  for (int row = 0; row < n_rows; row++) {
    for (int col = 0; col < n_cols; col++) {
      grid(row, col) = R::runif(0, 1);
    }
  }

  for (int it = 0; it < iterations; it++) {
    source_row = rand() % n_rows;
    source_col = rand() % n_cols;
    source_val = grid(source_row, source_col);
    span_row = rand() % max_span;
    span_col = rand() % max_span;
    row = source_row - span_row;
    col = source_col - span_col;
    do {
      c = wrap(col, n_cols);
      do {
        r = wrap(row, n_rows);
        grid(r, c) = (grid(r, c) + source_val) / 2;
        row++;
      } while (row < source_row + span_row);
      col++;
    } while (col < source_col + span_col);
  }

  return grid;
}
