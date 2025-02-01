defmodule Wavelets.DWTTest do
  use ExUnit.Case
  doctest Wavelets.DWT

  import Wavelets.TestHelpers
  alias Wavelets.{DWT, Filters}

  test "correctly transforms signal with Haar wavelet" do
    signal = [4.0, 4.0, 2.0, 2.0]
    # Haar wavelet
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    # Expected values for Haar wavelet
    expected_approx = [4.0 * :math.sqrt(2), 2.0 * :math.sqrt(2)]
    expected_details = [0.0, 0.0]

    assert_close(approx, expected_approx)
    assert_close(details, expected_details)
  end

  test "preserves energy" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    filter = Filters.Daubechies.get(2)

    {approx, details} = DWT.forward_1d(signal, filter)

    original_energy = signal |> Enum.map(&(&1 * &1)) |> Enum.sum()
    transform_energy = (approx ++ details) |> Enum.map(&(&1 * &1)) |> Enum.sum()

    assert_in_delta original_energy, transform_energy, 1.0e-10
  end

  test "preserves energy in 2D transform" do
    signal = [
      [1.0, 2.0, 3.0, 4.0],
      [5.0, 6.0, 7.0, 8.0],
      [9.0, 10.0, 11.0, 12.0],
      [13.0, 14.0, 15.0, 16.0]
    ]

    filter = Filters.Daubechies.get(1)

    {approx, {h_details, v_details, d_details}} = DWT.forward_2d(signal, filter)

    original_energy =
      signal |> List.flatten() |> Enum.map(&(&1 * &1)) |> Enum.sum()

    transform_energy =
      (List.flatten(approx) ++
         List.flatten(h_details) ++
         List.flatten(v_details) ++ List.flatten(d_details))
      |> Enum.map(&(&1 * &1))
      |> Enum.sum()

    assert_in_delta original_energy, transform_energy, 1.0e-10
  end
end
