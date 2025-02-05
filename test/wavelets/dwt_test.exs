defmodule Wavelets.DWTTest do
  use ExUnit.Case
  doctest Wavelets.DWT

  import Wavelets.TestHelpers
  alias Wavelets.{DWT, Filters}

  test "correctly transforms signal with Haar wavelet" do
    signal = [4.0, 4.0, 2.0, 2.0]
    filter = Filters.Daubechies.get(1)
    {approx, details} = DWT.forward_1d(signal, filter)

    # Expected values matching PyWavelets orthonormal scaling
    expected_approx = [4.0 * :math.sqrt(2), 2.0 * :math.sqrt(2)]
    expected_details = [0.0, 0.0]

    assert_close(approx, expected_approx)
    assert_close(details, expected_details)
  end

  test "preserves energy" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    filter = Filters.Daubechies.get(2)

    {approx, details} = DWT.forward_1d(signal, filter)

    original_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    transform_energy =
      Enum.sum(Enum.map(approx, &(&1 * &1))) +
        Enum.sum(Enum.map(details, &(&1 * &1)))

    # Energy should be exactly preserved
    assert_in_delta original_energy, transform_energy, 1.0e-8
  end
end
