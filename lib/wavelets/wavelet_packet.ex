defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction.
  """

  alias Wavelets.{DWT, IDWT, Filter}

  @doc """
  Performs 1D wavelet packet decomposition.
  Creates a tree structure where each node contains coefficients at different
  scales and frequency bands.
  """
  @spec decompose_1d(
          list(number),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) ::
          map()
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    tree = %{
      {0, 0} => signal
    }

    do_decompose_1d(tree, filter, levels)
  end

  @doc """
  Reconstructs the signal from its wavelet packet decomposition.
  """
  @spec reconstruct_1d(map(), Filter.t(), Wavelets.precision()) :: list(number)
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    do_reconstruct_1d(tree, filter, max_level)
  end

  # Private helper functions
  defp do_decompose_1d(tree, _filter, 0), do: tree

  defp do_decompose_1d(tree, filter, level) do
    new_nodes =
      for j <- 0..(:math.pow(2, level - 1) |> trunc() |> Kernel.-(1)),
          node = {level - 1, j},
          Map.has_key?(tree, node) do
        signal = Map.get(tree, node)
        {approx, details} = DWT.forward_1d(signal, filter)

        %{
          {level, j * 2} => approx,
          {level, j * 2 + 1} => details
        }
      end

    new_tree =
      Enum.reduce(new_nodes, tree, fn node_map, acc ->
        Map.merge(acc, node_map)
      end)

    do_decompose_1d(new_tree, filter, level - 1)
  end

  defp do_reconstruct_1d(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end

  defp do_reconstruct_1d(tree, filter, level) do
    new_nodes =
      for j <- 0..(:math.pow(2, level - 2) |> trunc() |> Kernel.-(1)) do
        approx = Map.get(tree, {level, j * 2})
        details = Map.get(tree, {level, j * 2 + 1})
        reconstructed = IDWT.inverse_1d(approx, details, filter)

        {{level - 1, j}, reconstructed}
      end

    new_tree =
      Enum.reduce(new_nodes, tree, fn {node, signal}, acc ->
        Map.put(acc, node, signal)
      end)

    do_reconstruct_1d(new_tree, filter, level - 1)
  end
end
