ExUnit.start()

defmodule Wavelets.TestHelpers do
  import ExUnit.Assertions

  def assert_close(list1, list2, tolerance \\ 1.0e-10) do
    assert length(list1) == length(list2)

    Enum.zip(list1, list2)
    |> Enum.each(fn {a, b} -> assert_in_delta(a, b, tolerance) end)
  end

  def assert_close_2d(matrix1, matrix2, tolerance \\ 1.0e-10) do
    assert length(matrix1) == length(matrix2)
    assert length(hd(matrix1)) == length(hd(matrix2))

    Enum.zip(matrix1, matrix2)
    |> Enum.each(fn {row1, row2} -> assert_close(row1, row2, tolerance) end)
  end
end
