defmodule Wavelets.FilterTest do
  use ExUnit.Case
  doctest Wavelets.Filter

  alias Wavelets.{Filter, Filters}

  describe "Daubechies filters" do
    test "creates valid Haar filter" do
      filter = Filters.Daubechies.get(1)

      assert filter.name == "db1"
      assert filter.family == :daubechies
      assert filter.vanishing_moments == 1
      assert length(filter.decomposition_low_pass) == 2
    end

    test "creates valid Daubechies-4 filter" do
      filter = Filters.Daubechies.get(2)

      assert filter.name == "db2"
      assert filter.family == :daubechies
      assert filter.vanishing_moments == 2
      assert length(filter.decomposition_low_pass) == 4
    end
  end

  describe "Symlet filters" do
    test "creates valid Symlet-4 filter" do
      filter = Filters.Symlet.get(4)

      assert filter.name == "sym4"
      assert filter.family == :symlet
      assert filter.vanishing_moments == 4
      assert length(filter.decomposition_low_pass) == 8
    end

    test "creates valid Symlet-8 filter" do
      filter = Filters.Symlet.get(8)

      assert filter.name == "sym8"
      assert filter.family == :symlet
      assert filter.vanishing_moments == 8
      assert length(filter.decomposition_low_pass) == 16
    end
  end

  describe "Morlet filter" do
    test "creates valid Morlet filter with default parameters" do
      filter = Filters.Morlet.get()

      assert filter.name == "morl"
      assert filter.family == :morlet
      assert is_nil(filter.vanishing_moments)
      assert length(filter.decomposition_low_pass) > 0
    end

    test "creates valid Morlet filter with custom sigma" do
      sigma = 3.0
      filter = Filters.Morlet.get(sigma)

      assert filter.support_width == trunc(8 * sigma)

      assert length(filter.decomposition_low_pass) ==
               2 * filter.support_width + 1
    end
  end
end
