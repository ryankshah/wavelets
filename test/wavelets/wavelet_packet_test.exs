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

    # Verify root node contains original signal
    assert_close(Map.get(tree, {0, 0}), signal)
  end

  test "perfectly reconstructs signal" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    filter = Filters.Daubechies.get(2)
    levels = 3

    tree = WaveletPacket.decompose_1d(signal, filter, levels)
    reconstructed = WaveletPacket.reconstruct_1d(tree, filter)

    # Use a slightly larger tolerance for numerical stability
    assert_close(signal, reconstructed, 1.0e-8)
  end

  test "handles different wavelet filters" do
    signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    filters = [
      Filters.Daubechies.get(1),
      Filters.Daubechies.get(2),
      Filters.Symlet.get(4)
    ]

    Enum.each(filters, fn filter ->
      tree = WaveletPacket.decompose_1d(signal, filter, 2)
      reconstructed = WaveletPacket.reconstruct_1d(tree, filter)
      assert_close(signal, reconstructed, 1.0e-8)

      # Also verify energy conservation
      original_energy = Enum.sum(Enum.map(signal, &(&1 * &1)))

      tree_energy =
        tree
        |> Map.values()
        |> Enum.map(&Enum.sum(Enum.map(&1, fn x -> x * x end)))
        |> Enum.sum()
        # Average over all nodes
        |> Kernel./(Map.size(tree))

      assert_in_delta original_energy, tree_energy, 1.0e-8
    end)
  end
end
