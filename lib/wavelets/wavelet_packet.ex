defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction using base DWT
  """

  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition
  """
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    # Each node stores normalized coefficients
    tree = %{{0, 0} => signal}

    # Build each level using normalized DWT
    Enum.reduce(1..levels, tree, fn level, acc ->
      build_level(acc, filter, level)
    end)
  end

  defp build_level(tree, filter, level) do
    prev_level = level - 1

    # Get nodes from previous level
    prev_nodes =
      for {key = {^prev_level, _}, signal} <- tree, into: %{} do
        {key, signal}
      end

    # Create new nodes with DWT
    new_nodes =
      Enum.flat_map(Map.to_list(prev_nodes), fn {{_, j}, signal} ->
        {approx, details} = DWT.forward_1d(signal, filter)

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
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    reconstruct_level(tree, filter, max_level)
  end

  defp reconstruct_level(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end

  defp reconstruct_level(tree, filter, level) do
    # Reconstruct pairs of nodes at current level
    nodes =
      for j <- 0..(:math.pow(2, level - 1) |> trunc() |> Kernel.-(1)),
          approx = Map.get(tree, {level, j * 2}),
          details = Map.get(tree, {level, j * 2 + 1}),
          approx != nil and details != nil do
        {{level - 1, j}, IDWT.inverse_1d(approx, details, filter)}
      end
      |> Map.new()

    # Merge with existing tree and continue reconstruction
    Map.merge(tree, nodes)
    |> reconstruct_level(filter, level - 1)
  end
end
