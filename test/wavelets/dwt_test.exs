defmodule Wavelets.DWTTest do
  use ExUnit.Case
  doctest Wavelets.DWT

  import Wavelets.TestHelpers
  alias Wavelets.{DWT, Filter, Filters}

  describe "forward_1d/3" do
    test "correctly transforms signal with Haar wavelet" do
      signal = [4.0, 4.0, 2.0, 2.0]
      # Haar wavelet
      filter = Filters.Daubechies.get(1)

      {approx, details} = DWT.forward_1d(signal, filter)

      expected_approx = [5.656854249492381, 2.82842712474619]
      expected_details = [0.0, 0.0]

      assert_close(approx, expected_approx)
      assert_close(details, expected_details)
    end

    test "preserves energy" do
      signal = generate_test_signal()
      filter = Filters.Daubechies.get(2)

      {approx, details} = DWT.forward_1d(signal, filter)

      original_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))
      transform_energy = Enum.sum(Enum.map(approx ++ details, &(&1 * &1)))

      assert_in_delta original_energy, transform_energy, 1.0e-10
    end
  end

  describe "forward_2d/3" do
    test "correctly decomposes 2D signal" do
      signal = generate_test_signal_2d()
      filter = Filters.Daubechies.get(1)

      {approx, {h_details, v_details, d_details}} =
        DWT.forward_2d(signal, filter)

      assert length(approx) == div(length(signal), 2)
      assert length(h_details) == div(length(signal), 2)
      assert length(v_details) == div(length(signal), 2)
      assert length(d_details) == div(length(signal), 2)
    end

    test "preserves energy in 2D transform" do
      signal = generate_test_signal_2d()
      filter = Filters.Daubechies.get(2)

      {approx, {h_details, v_details, d_details}} =
        DWT.forward_2d(signal, filter)

      original_energy =
        signal
        |> List.flatten()
        |> Enum.map(&(&1 * &1))
        |> Enum.sum()

      transform_energy =
        (List.flatten(approx) ++
           List.flatten(h_details) ++
           List.flatten(v_details) ++
           List.flatten(d_details))
        |> Enum.map(&(&1 * &1))
        |> Enum.sum()

      assert_in_delta original_energy, transform_energy, 1.0e-10
    end
  end
end
