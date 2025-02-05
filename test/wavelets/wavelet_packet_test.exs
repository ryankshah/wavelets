defmodule Wavelets.WaveletPacketTest do
  use ExUnit.Case
  doctest Wavelets.WaveletPacket

  import Wavelets.TestHelpers
  alias Wavelets.{Filters, WaveletPacket}

  test "creates correct tree structure" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    filter = Filters.Daubechies.get(1)
    levels = 2

    tree = WaveletPacket.decompose_1d(signal, filter, levels)

    # Check tree structure
    assert Map.has_key?(tree, {0, 0})
    assert Map.has_key?(tree, {1, 0})
    assert Map.has_key?(tree, {1, 1})
    assert Map.has_key?(tree, {2, 0})
    assert Map.has_key?(tree, {2, 1})
    assert Map.has_key?(tree, {2, 2})
    assert Map.has_key?(tree, {2, 3})

    # Root node should contain original signal
    assert Map.get(tree, {0, 0}) == signal
  end

  test "perfectly reconstructs signal" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    filter = Filters.Daubechies.get(2)
    levels = 3

    tree = WaveletPacket.decompose_1d(signal, filter, levels)
    reconstructed = WaveletPacket.reconstruct_1d(tree, filter)

    assert_close(signal, reconstructed, 1.0e-8)
  end

  test "handles different wavelet filters" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    original_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

    filters = [
      Filters.Daubechies.get(1),
      Filters.Daubechies.get(2),
      Filters.Symlet.get(4)
    ]

    Enum.each(filters, fn filter ->
      tree = WaveletPacket.decompose_1d(signal, filter, 2)
      reconstructed = WaveletPacket.reconstruct_1d(tree, filter)
      assert_close(signal, reconstructed, 1.0e-8)

      # Check energy conservation at each level
      Enum.each(0..2, fn level ->
        level_coeffs = for j <- 0..(trunc(:math.pow(2, level)) - 1),
                          coeffs = Map.get(tree, {level, j}),
                          coeffs != nil,
                          do: coeffs |> List.flatten()
        
        level_energy = level_coeffs |> List.flatten() |> Enum.map(&(&1 * &1)) |> Enum.sum()
        assert_in_delta original_energy, level_energy, 1.0e-8
      end)
    end)
  end
end