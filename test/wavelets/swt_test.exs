defmodule Wavelets.SWTTest do
  use ExUnit.Case
  doctest Wavelets.SWT

  import Wavelets.TestHelpers
  alias Wavelets.Filters
  alias Wavelets.SWT

  describe "transform_1d/4" do
    test "computes correct number of levels" do
      signal = generate_test_signal()
      filter = Filters.Daubechies.get(1)
      levels = 3

      {approximations, details} = SWT.transform_1d(signal, filter, levels)

      # Including original signal
      assert length(approximations) == levels + 1
      assert length(details) == levels
    end

    test "preserves signal length at each level" do
      signal = generate_test_signal()
      filter = Filters.Daubechies.get(1)
      levels = 2

      {approximations, details} = SWT.transform_1d(signal, filter, levels)

      signal_length = length(signal)

      Enum.each(approximations, fn approx ->
        assert length(approx) == signal_length
      end)

      Enum.each(details, fn detail ->
        assert length(detail) == signal_length
      end)
    end

    test "handles different wavelet filters" do
      signal = generate_test_signal()

      filters = [
        Filters.Daubechies.get(1),
        Filters.Daubechies.get(2),
        Filters.Symlet.get(4)
      ]

      Enum.each(filters, fn filter ->
        {approximations, details} = SWT.transform_1d(signal, filter, 2)
        assert length(approximations) == 3
        assert length(details) == 2
      end)
    end
  end

  # Helper function for generating test signals
  defp generate_test_signal do
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  end
end
