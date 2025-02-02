defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction
  that maintains orthonormality and energy preservation
  """

  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition ensuring orthonormality
  """
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    # Start with root node containing original signal
    base_tree = %{{0, 0} => signal}

    # Build each level
    Enum.reduce(1..levels, base_tree, fn level, tree ->
      build_level(tree, filter, level)
    end)
  end

  defp build_level(tree, filter, level) do
    prev_level = level - 1

    # Get nodes from previous level
    prev_nodes =
      for {key = {^prev_level, _}, signal} <- tree, into: %{} do
        {key, signal}
      end

    # Create new nodes using DWT (which now has correct scaling)
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
  ensuring perfect reconstruction
  """
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    # Reconstruct from the highest level down
    reconstruct_level(tree, filter, max_level)
  end

  defp reconstruct_level(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end

  defp reconstruct_level(tree, filter, level) do
    # Process all node pairs at current level
    nodes =
      for j <- 0..(:math.pow(2, level - 1) |> trunc() |> Kernel.-(1)),
          approx = Map.get(tree, {level, j * 2}),
          details = Map.get(tree, {level, j * 2 + 1}),
          approx != nil and details != nil do
        # IDWT will handle proper scaling
        {{level - 1, j}, IDWT.inverse_1d(approx, details, filter)}
      end
      |> Map.new()

    # Continue reconstruction at next level up
    Map.merge(tree, nodes)
    |> reconstruct_level(filter, level - 1)
  end
end
