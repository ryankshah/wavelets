defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction with proper normalization
  """

  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition
  """
  def decompose_1d(signal, filter, levels, precision \\ :double) do
    # Initialize tree with original signal
    tree = %{{0, 0} => signal}

    # Build tree with proper orthonormal scaling at each level
    Enum.reduce(1..levels, tree, fn level, acc ->
      build_level(acc, filter, level, precision)
    end)
  end

  defp build_level(tree, filter, level, precision) do
    prev_level = level - 1

    # Get nodes from previous level
    prev_nodes =
      for {key = {^prev_level, _}, signal} <- tree, into: %{} do
        {key, signal}
      end

    # Create new nodes with proper scaling
    new_nodes =
      prev_nodes
      |> Enum.flat_map(fn {{_, j}, signal} ->
        {approx, details} = DWT.forward_1d(signal, filter, precision)

        [
          {{level, j * 2}, approx},
          {{level, j * 2 + 1}, details}
        ]
      end)
      |> Map.new()

    Map.merge(tree, new_nodes)
  end

  @doc """
  Reconstructs signal from wavelet packet decomposition
  """
  def reconstruct_1d(tree, filter, precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    reconstruct_level(tree, filter, max_level, precision)
  end

  defp reconstruct_level(tree, _filter, 0, _precision) do
    Map.get(tree, {0, 0})
  end

  defp reconstruct_level(tree, filter, level, precision) do
    nodes =
      for j <- 0..(:math.pow(2, level - 1) |> trunc() |> Kernel.-(1)),
          approx = Map.get(tree, {level, j * 2}),
          details = Map.get(tree, {level, j * 2 + 1}),
          approx != nil and details != nil do
        reconstructed = IDWT.inverse_1d(approx, details, filter, precision)
        {{level - 1, j}, reconstructed}
      end
      |> Map.new()

    Map.merge(tree, nodes)
    |> reconstruct_level(filter, level - 1, precision)
  end
end
