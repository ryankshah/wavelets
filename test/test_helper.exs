ExUnit.start()

defmodule Wavelets.TestHelpers do
  import ExUnit.Assertions

  def assert_close(list1, list2, tolerance \\ 1.0e-10) do
    assert length(list1) == length(list2)
    
    Enum.zip(list1, list2)
    |> Enum.each(fn {a, b} ->
      assert_in_delta(a, b, tolerance)
    end)
  end

  def assert_close_2d(matrix1, matrix2, tolerance \\ 1.0e-10) do
    assert length(matrix1) == length(matrix2)
    assert length(hd(matrix1)) == length(hd(matrix2))
    
    Enum.zip(matrix1, matrix2)
    |> Enum.each(fn {row1, row2} ->
      assert_close(row1, row2, tolerance)
    end)
  end

  def generate_test_signal(length \\ 8) do
    1..length |> Enum.map(&(&1 * 1.0))
  end

  def generate_test_signal_2d(size \\ 4) do
    for i <- 1..size do
      for j <- 1..size do
        (i - 1) * size + j * 1.0
      end
    end
  end
end